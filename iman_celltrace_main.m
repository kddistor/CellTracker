%IMAN_CELLTRACE_MAIN

function [r, ri, GMD, xtr] = iman_celltrace_main(p, r, ri)    

%Check for existing reader (e.g. from previous or failed run)
use_old_r = exist('r','var')  && ~isempty(r)  && ...
            exist('ri','var') && ~isempty(ri);

%Extract parameters for parallel runs from structure
%FIXME include defaults for robustness in calling
nW = p.nW;  xtalk = p.xtalk;  bval = p.bval;  bk = p.bk;  op = p.op;
dispsamp = p.disp;
thet = 0:(2*pi/119):(2*pi);

%Determine segmentation procedure
switch regexprep(op.local, ...
        {'(cyt.*)|(^c$)', '(nuc.*)|(^n$)'},{'cyto','nucl'})
    case 'cyto' %Cytoplasmic signals
        op.iscyto = true;
    case 'nucl' %Nuclear signals
        op.iscyto = false;
end
        
%Parallelization setup
if nW > 1 && matlabpool('size') == 0; matlabpool('local', nW); end

runid = ceil(rand*1e4);                     %Assign random run ID number
tsn = [p.tname,'\IMAN_TEMP_',num2str(runid),'_n'];%Temporary Save Name


%% DATA PREPARATION  -----------------------------------------------------
%Establish Reader for ND2 file
if ~use_old_r
    msg = sprintf(['Loading ND2 file Reader.  ',...
        'Start time is:  ', datestr(now, 'mm/dd HH:MM PM')]);    disp(msg);
    [r, ri] = bfopen_custom(p.fname);   
    msg = sprintf(['Reader loaded.  ', datestr(now, 'mm/dd HH:MM PM')]);
    disp(msg);
else
    msg = sprintf('Using existing ND2 file Reader (r, ri) provided.');
    disp(msg);
end
msg = [];   %Clear msg to prevent partial erasing on next update

%Extract MetaData from ND2 file
%   Will warn if missing metadata and prompt for input as needed.
GMD = nd2_meta(r, ri, p.bkmd);
%   Double check GMD for flaws (e.g. out of range values?)

%Enforce backup values as requested
if ~isempty(p.bkenforce);
    for s = 1:size(p.bkenforce,1)
        GMD.(p.bkenforce{s,1}).(p.bkenforce{s,2}) = ...
            p.bkmd.(p.bkenforce{s,1}).(p.bkenforce{s,2}); 
    end
end

%FIX CORRUPTED ND2 FILE ORDER?  (via .lblr)
%   Create pseudo-Series per XY position?

%Determine XY storage order
if numel(ri) >= max(p.xypos) && ri(1).ZTCsize(1) == 1
    xySeries = true;    sz = 1;        
elseif ri(1).ZTCsize(1) > max(p.xypos)
    xySeries = false;   r.setSeries(0);
else
    %FIXME Some sort of error I suppose
    warning(['XY storage order could not be determined. Aborting', ...
             ' and returning ND2 Reader and MetaData']);
    xtr = [];  return;
end

%Sample display parameters
dsx = min(640,GMD.cam.PixNumX);     dsy = min(540,GMD.cam.PixNumY);  
dscx = dsx./GMD.cam.PixNumX;        dscy = dsy./GMD.cam.PixNumY;


%% PREPROCESSING PREPARATION  --------------------------------------------
%Channel ID definition (for use below)
nChan = r.getSizeC;
cid = cell2struct(num2cell(1:nChan), GMD.exp.Channel, 2);
%   Ensure channel reference uniformity
cn = fieldnames(cid);
%   Prepare cross talk correction if applicable
xtalkcorr = false;  xtr = [];
if ~isempty(xtalk.D)
    xtalk = structfun(@(x) cn{ ~cellfun( 'isempty', regexpi(cn, ...
        x, 'start')) }, xtalk, 'UniformOutput', false);
    %   Get cross-talk ratios for simple removal
    xtr = iman_xtalk_correction(xtalk.I, xtalk.D, GMD, []);
    xtalkcorr = true;
end

%Use zero (0) as placeholder for 1 digit xypos
dhold = '0';    

%Establish memory constraints
[uV,sV] = memory;  maxmem = sV.PhysicalMemory.Available;  %#ok<ASGLU> %Available memory

%Sequentially process each xy position desired
t_tot = tic;  nxy = numel(p.xypos);  t_perxy = [];  emsg = [];  errid = 0;
for h = 1:nxy
    t_xy = tic;
    try  %Protect other xy runs from failures
    %Ensure memory is not occupied from previous run
    clear imstack valcube masks movieInfo tracksFinal coord
        
    %Select current Series in Reader, or Z slice (to index xy)
    if xySeries; r.setSeries(p.xypos(h)-1);  else   sz = p.xypos(h);   end
    nT = r.getSizeT;        %Get number of time points for this XY
    
    %Get background intensities, unless dynamic
    %FIXME Implement 'global' background from empty well
    if bk(h).fix    %IF fixed and provided
        bki = bk(h).reg - bval;
    elseif ~bk(h).dyn   %IF static and region provided for image 1
        im = zeros(r.getSizeY, r.getSizeX, nChan);
        for sc = 1:nChan
            im(:,:,sc) = bfGetPlane( r, ri(p.xypos(h)).lblr(sz,1,sc) );
        end
        %Refine image for accurate background estimate
        im = iman_baseline_removal(im, bval);
        %Correct objective view bias (collection efficiency)
        [im, e_inv] = iman_objective_correction(im, GMD); %#ok<NASGU>
        %Correct cross-talk (optional)
        if xtalkcorr
            %   Perform cross-talk subtraction (here, stacked channels)
            im(:,:,cid.(xtalk.D)) = im(:,:,cid.(xtalk.D)) ...
                - xtr.(xtalk.D).(xtalk.I).*im(:,:,cid.(xtalk.I));
        end
        bki = mean(mean( ...
            im( bk(h).reg(3):bk(h).reg(4), ...
            bk(h).reg(1):bk(h).reg(2), : ), 1 ), 2);
    else bki = [];
    end
    
    %Check memory requirement to load nW series (4 Bytes for single)
    MemLoad = nT*nChan*r.getSizeX*r.getSizeY*4*10;  %Safety factor of 10
    %   Determine how to slice
    nP  = ceil( MemLoad/maxmem );               %Number of Parts
    ipw = ceil( (nT)./(nP.*nW) );               %Images per Worker
    %Define slicing for parallel processing (per part)
    nper = ipw*ones(1,nP*nW);  dif = nP*nW*ipw - nT - 1; %Fix overestimate
    nper(end-dif:end) = nper(end-dif:end) - 1;  %Minus 1 excess per Worker
    
    %% Process images in parallel batches
    mTEMP = cell(1,nP);
    for sn = 1:nP
    %Load Images
        %Display status message
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
        [msg, emsg] = update_timer(t_el, t_rem, msg, emsg);
        msg = sprintf(['Loading data for part %d/%d of XY ', ...
            '%d (%d/%d).'], sn, nP, p.xypos(h), h, nxy);    disp(msg);
    %Fill data for this Batch
    imstack = cell(1,nW);
    for sw = 1:nW   %Fill Worker
        wnum = (sn-1)*nW + sw;  %Worker number (from start)
        for sp = 1:nper(wnum)   %Fill Part (indexing into Time)
            st = sum(nper(1:wnum-1)) + sp;  %Get Time index
            for sc = 1:nChan                %Fill Channels
                imstack{sw}{sp}(:,:,sc) = ...
                    bfGetPlane( r, ri(p.xypos(h)).lblr(sz,st,sc) );
            end
        end
    end      
    
    %Process images in parallel (each time point)
        %Display status message
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
        [msg, emsg] = update_timer(t_el, t_rem, msg, emsg);
        msg = sprintf(['Identifying cells for part %d/%d of XY ', ...
            '%d (%d/%d).'], sn, nP, p.xypos(h), h, nxy);    disp(msg);
    movieInfo = cell(1,nW);  bkis = cell(nW,1);
    parfor ps = 1:nW   %Parallel Loop
        e_inv = [];  %Initialize temporary variable
        for s = 1:length(imstack{ps})
            %Refine Image
            im = double(imstack{ps}{s});   %Get current image (as a double)
            %Remove baseline value from image
            %   The scalar addition should be removed prior to corrections
            im = iman_baseline_removal(im, bval);
            %Correct objective view bias (collection efficiency)
            if s == 1; [im, e_inv] = iman_objective_correction(im, GMD);
            else    im = iman_objective_correction(im, GMD, e_inv);   end
            %Correct cross-talk (optional)
            if xtalkcorr
            %   Perform cross-talk subtraction (here, stacked channels)
            im(:,:,cid.(xtalk.D)) = im(:,:,cid.(xtalk.D)) ...
                - xtr.(xtalk.D).(xtalk.I).*im(:,:,cid.(xtalk.I)); %#ok<PFBNS>
            end
            
            if dispsamp && ps==1 && s==1;   %Store for display only
                im_or{ps} = imresize(imstack{ps}{s}, [dsy,dsx]/2); end
            %Store refined imaged
            imstack{ps}{s} = single(im);
            %Get background image values
            if bk(h).dyn && ~bk(h).fix   %#ok<PFBNS>
                bkit = mean(mean( ...
                    im( bk(h).reg(3):bk(h).reg(4), ...
                        bk(h).reg(1):bk(h).reg(2), : ), 1 ), 2);
                bkis{ps}(s,:) = bkit;    %Store background used
            else bkit = bki;    
            end
            %Perform Cell Idenification
            movieInfo{ps}{1,s} = iman_cellid(im, op, bkit);
        end 
    end
        %Periodic sample display -----------------------------------------
        if dispsamp
            ps = 1; s = 1;
            %SHOW REFINED IMAGES
            figure(runid); clf reset; %Figure with unique ID for this run
            set(runid, 'Position', [80, 50, 700, 920]);
            for ds = 1:nChan
                subplot(nChan,2,1+(ds-1)*2);  
                set(gca,'Position', get(gca,'Position') + ...
                    [-0.08, -0.02*ds, 0.1, 0.06] );
                imshow(im_or{ps}(:,:,ds), [],...
                    'Border', 'tight', 'InitialMagnification', 'fit');
                lh = ylabel(GMD.exp.Channel{ds}); 
                lhp = get(lh,'Position'); ax = axis;
                set(lh, 'Position', [-0.025*ax(4), lhp(2:end)]);
                if ds == 1; title('Original'); end 
                subplot(nChan,2,2+(ds-1)*2); 
                set(gca,'Position', get(gca,'Position') + ...
                    [-0.05, -0.02*ds, 0.1, 0.06] );
                imshow(imresize(imstack{ps}{s}(:,:,ds), [dsy,dsx]/2), [], ...
                    'Border', 'tight', 'InitialMagnification', 'fit');
                if ds == 1; title('Refined'); end
            end
            %SHOW CELL ID RESULTS
            figure(runid+1); clf reset; %Figure with unique ID for this run
            set(runid+1, 'Position', [80, 50, 710, 600]);
            %Show Image
            imshow( imresize(imstack{ps}{s}(:,:,op.chan), [dsy,dsx]), [],...
                'Border', 'tight', 'InitialMagnification', 'fit');  hold on;
            %Plot centers
            scatter(movieInfo{ps}{1,s}.xCoord.*dscx, ...
                movieInfo{ps}{1,s}.yCoord.*dscy, 75, 'r+');
            %Plot circles representing ~nuclear area
            for ds = 1:length(movieInfo{ps}{1,s}.xCoord)
                ampy = sqrt(movieInfo{ps}{1,s}.amp(ds).*dscx.*dscy ./ pi);
                plot(movieInfo{ps}{1,s}.xCoord(ds).*dscx + ampy.*cos(thet), ...
                    movieInfo{ps}{1,s}.yCoord(ds).*dscy + ampy.*sin(thet), 'b:');
            end
            %Plot reference cell areas
            plot( (op.maxD/2 + 0.5*op.maxD.*cos(thet)).*dscx, ...
                (op.maxD/2 + 0.5*op.maxD.*sin(thet)).*dscy, 'r:' );
            plot( (op.maxD/2 + 0.5.*op.minD.*cos(thet)).*dscx, ...
                (op.maxD/2 + 0.5.*op.minD.*sin(thet)).*dscy, 'r:' );
            drawnow;
        end
        %Periodic sample display -----------------------------------------
    movieInfo = cat(2, movieInfo{:});     %Cat parallel layer
    mTEMP{1,sn} = cat(2, movieInfo{:});   %Cat serial layer
    %Save temp data locally for rapid reloading
    save([tsn,num2str(sn)], 'imstack');
    end
    movieInfo = cat(2, mTEMP{:});  clear mTEMP;
    
        %Update timing and display    
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
        [msg, emsg] = update_timer(t_el, t_rem, msg, emsg);
        msg = [];  emsg = [];  %Clear to prevent deleting scripttrack msg.
    
    
    %% Perform Cell Tracking
    
    %If a known shift in XY position occurred, provide magnitude and timing
    if  isfield(p,'xyshift') && ~isempty(p.xyshift)
        movieInfo = iman_xyshift(movieInfo, p.xyshift, -1);
    end
    
    %   For scriptTrackGeneral, movieInfo fields must include an ending
    %       column vector of zeros
    for s = 1:numel(movieInfo)
        movieInfo(s) = structfun(@(x)[x, zeros(size(x,1),1)], movieInfo(s),...
            'UniformOutput', false);
    end
    tracksFinal = scripttrackwrap(movieInfo, p);
    
    %Perform and track cleanup to deliver xC and yC matrices
    %   Provide sliced matrices for contouring
    coord = iman_trackcoords(tracksFinal);
    %If coords is not full length (no tracks are full), pad with NaN
    coord(:,(end+1):nT,:) = NaN;
    
    %Reverse any XY shift define
    if  isfield(p,'xyshift') && ~isempty(p.xyshift)
        coord = iman_xyshift(coord, p.xyshift, 1);
    end
    
    %Slice coordinate time to match imstacks
    coord = mat2cell(coord, size(coord,1), nper, 2);    %To nW*nP slices
    coord = mat2cell(coord, 1, nW*ones(nP,1));          %To nP slices
    
        %Update timing and display    
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
        disp(' ');  %Make up space for stripttrack printing
        [msg, emsg] = update_timer(t_el, t_rem, msg, emsg);
        

    %% Process contour definitions in parallel batches
    %Establish data output sequence (to be stored for reference)
    vcorder = iman_cellmask(double(imstack{1}{1}), [], [], op, []); %#ok<NASGU>
    %Load temp data from local storage
    vTEMP = cell(1,nP);  if bk(h).dyn && ~bk(h).fix; bTEMP = cell(nP,1); end
    mTEMP = cell(1,nP);
    for sn = 1:nP
        %Display status message
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h);
        [msg, emsg] = update_timer(t_el, t_rem, msg, emsg);
        msg = sprintf(['Masking and Contouring for part %d/%d of XY ', ...
            '%d (%d/%d).'], sn, nP, p.xypos(h), h, nxy);    disp(msg);
    %IF more than 1 part used, load part. Contains 'imstack' and 'tind'
    if nP > 1;    load([tsn,num2str(sn)]);    end 
    coord_t = coord{sn};
    valcube = cell(1,nW);
    masks = cell(1,nW);
    parfor ps = 1:nW   %Parallel Loop
        for s = 1:length(imstack{ps})
            %Perform Contour ID and Feature Extraction
            im = double(imstack{ps}{s});    %Get current image (as double)
            %Get background image values
            if bk(h).dyn && ~bk(h).fix;  bkit = bkis{ps}(s,:);   %#ok<PFBNS>
            else bkit = bki;  end
            [valcube{ps}{s}, masks{ps}{s}] = iman_cellmask(im, ...
                coord_t{ps}(:,s,:), tracksFinal, op, ...
                bkit);
        end
    end
        %Periodic sample display -------------------------------------
        if dispsamp 
            ps = 1; s = 1;
            figure(runid+2); clf reset; %Figure with unique ID for this run
            set(runid+2, 'Position', [80, 50, 710, 600]);
            %Shade image over nuclear and cytoplamsic masks
            mnv = min(min(imstack{ps}{s}(:,:,op.chan)));   
            mxv = max(max(imstack{ps}{s}(:,:,op.chan)));
            gsi = (imstack{ps}{s}(:,:,op.chan)-mnv)./(mxv-mnv);
            rsi = gsi;  bsi = rsi;
            for ds = 1:length(masks{ps}{s})
                rsi(masks{ps}{s}(ds).nuc) = rsi(masks{ps}{s}(ds).nuc)/4;
                bsi(masks{ps}{s}(ds).nuc) = bsi(masks{ps}{s}(ds).nuc)/4;
                gsi(masks{ps}{s}(ds).nuc) = gsi(masks{ps}{s}(ds).nuc)*3;
                rsi(masks{ps}{s}(ds).cyt) = rsi(masks{ps}{s}(ds).cyt)/4;
                gsi(masks{ps}{s}(ds).cyt) = gsi(masks{ps}{s}(ds).cyt)/4;
                bsi(masks{ps}{s}(ds).cyt) = bsi(masks{ps}{s}(ds).cyt)*3;
            end
            imshow( imresize(cat(3,rsi,gsi,bsi), [dsy,dsx]), [], ...
                'Border', 'tight', 'InitialMagnification', 'fit');
            clear rsi gsi bsi;  drawnow;
        end
        %Periodic sample display -------------------------------------
    valcube = cat(2, valcube{:});   vTEMP{sn} = cat(2, valcube{:}); 
    masks = cat(2, masks{:});       mTEMP{sn} = cat(2, masks{:});
    if bk(h).dyn;   bTEMP{sn} = cat(1, bkis{:});    end
    
    end
    valcube = cat(2, vTEMP{:});     clear vTEMP;  %#ok<NASGU>
    masks = cat(2, mTEMP{:});       clear mTEMP;
    if bk(h).dyn;  bki = cat(1, bTEMP{:});  clear bTEMP;   end %#ok<NASGU>
    %Convert masks to single to save space
    masks = changeprecision(masks, 'single'); %#ok<NASGU>
    
    %May use an index to convert masks to a 'linear' structure with no
    %   empty cells.  masks(index(track,time));
    %   However, little savings are achieved for dense data.
    %       idex = nan(size(masks));
    %       goodi = ~cellfun('isempty',{masks.nuc});
    %       idex(goodi) = 1:nnz(goodi);
    %       masks = masks(goodi);
    %       Usage:  masks(idex(track, time));
    
    
    %% Finalize and Save    
    %Save Processed Data to desired located
    xystring = ['_xy', dhold(1:double(p.xypos(h)<10)), num2str(p.xypos(h))];
    save([p.sname, xystring, '.mat'], 'valcube', 'vcorder', 'bki', 'masks');
    %Clear temporary files from disk
    for sn = 1:nP;  delete([tsn,num2str(sn),'.mat']);  end
    catch ME  %ERROR CATCH OPERATIONS
        errid = errid + 1;
        %Display error message
        dmsg = sprintf(repmat('\b', 1, numel(msg)+numel(emsg)+2));  
        disp(dmsg);        
        msg = sprintf(['Processing FAILED on XY ', ...
            '%d (%d/%d).'], p.xypos(h), h, nxy);    disp(msg);
        disp(ME.message);  emsg = []; msg = [];  
        %Produce an Error Report and skip to the next XY position
        MElog{errid}(1:2) = { msg, [ME.identifier, ' - ', ME.message] }; %#ok<AGROW>
        for se = 1:length(ME.stack)
            MElog{errid}(2 + se) =   {sprintf('%s  - Line %d', ...
                ME.stack(se).file, ME.stack(se).line) };  %#ok<AGROW>
        end
        %End of MElog is a structure with relevant info
        MElog{errid}{end} = struct('h',h); %#ok<AGROW>
        continue
    end
    
    %Update timing and display
    t_el = toc(t_tot)/60;   t_perxy = t_el/h;   t_rem = t_perxy*(nxy-h);
    [msg, emsg] = update_timer(t_el, t_rem, msg, emsg);
end
%Save MetaData and XTalk information
save([p.sname, '_Global.mat'], 'GMD', 'xtr', 'p');

%Close the parallel enviroment and return
if matlabpool('size') > 0;    matlabpool close;   end

if errid > 0
    %Open Error log file
    ct = now;
    erfname = ['CellTrace_ErrorLog_',datestr(ct,30),'.txt'];
    erfid = fopen(erfname, 'w+');
    fprintf(erfid, '%s\n%s\n\n', ['Error Log File for iman_celltrace_main,',...
        ' run ', datestr(ct, 31), ' for file: '], p.fname);
    for se = 1:errid
        fprintf(erfid, '\nError %d\n', se);         %New title line
        fprintf(erfid, '\nIn XY Position %d\n', ...
            p.xypos(MElog{se}{end}.h));          %XY Pos
        for sse = 1 : (length(MElog{se}) - 1)    %Print error messages
            fprintf(erfid, '%s\n', MElog{se}{sse});
        end
    end
    fclose(erfid);
end

end


%Timer display update function
function [msg, emsg] = update_timer(t_el, t_rem, msg, emsg)
fprintf(repmat('\b', 1, numel(msg)  + ~isempty(msg) + ...
                        numel(emsg) + ~isempty(emsg)  ));
msg = [];
%IF t_rem not yet defined, skip time estimation
if isempty(t_rem);      return;     end
emsg = sprintf(['Time elapsed: %dh %dm.  Estimate time remaining: ',...
    '%dh %dm.'], floor(t_el/60),  round(mod(t_el,60)), ...
                 floor(t_rem/60), round(mod(t_rem,60))  );
disp(emsg);
end