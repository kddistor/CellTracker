%IMAN_CELLTRACE_PARBF
%   Main image analysis function for extracting cell trace data from ND2
%   file images, using a parallelized Bioformats reader to access image
%   data.

function [r, ri, GMD, xtr] = iman_celltrace_parbfccount(p, timeRange, varargin)

%% INPUT PARSING AND CHECKING  -------------------------------------------

%Default parameters
bk = struct('reg',[],'dyn',false,'fix',false,'altxy',[],'bkonly',false);
op = struct('chan',[],'local',[],'cname',[],'maxD',[],'minD',[],...
    'minF',[],'diag',false,'stormask',false);
xtalk = struct('I','YFP','D','RFP');

%Extract parameters for parallel runs from input structure (avoids copying)
fname = p.fname;    xypos = p.xypos;    nW = p.nW;
bval = p.bval;      dispsamp = p.disp;
nbk = numel(p.bk);  bk(2:nbk) = bk;
for s = fieldnames(p.bk)'; [bk.(s{1})] = deal(p.bk.(s{1})); end
for s = fieldnames(p.op)'; op.(s{1}) = p.op.(s{1}); end
for s = fieldnames(p.xtalk)'; xtalk.(s{1}) = p.xtalk.(s{1}); end

thet = 0:(2*pi/119):(2*pi);  %Used to plot display approx. nuclear location
nxy = numel(xypos);  %Number of XY positions to be run
ncname = numel(op.cname);

%Check on background information provided

if nbk > nxy; bk = bk(1:nxy);
    warning(['Extra background information found: numel(bk) > numel',...
        '(xypos). Removing last ',num2str(nbk - nxy),' element(s) from bk.']);
elseif nbk < nxy; error(['Fewer backgrounds provided than XY positions.',...
        ' The background structure (bk) must contain an element for ',...
        'each XY position to be used.']);
end
%Set indices for re-ordering run when some XYs depend on others' background
use_altXY = ~cellfun('isempty',{bk.altxy});  xyrev(xypos) = 1:nxy;
altXYi = NaN(1,nxy); altXYi(use_altXY) = xyrev([bk(use_altXY).altxy]);
%Check validity of re-ordering structure
bkvalid = use_altXY;
bkvalid(~use_altXY) = arrayfun(@(x)( x.fix && numel(x.reg) == ncname ) || ...
    ( ~x.fix && isequal(size(x.reg),[2,2]) ),   bk(~use_altXY));
unaltXYi = unique(altXYi(use_altXY));    %Check for mislabeled alts
bkvalid( unaltXYi(use_altXY(unaltXYi)) ) = false;
if ~all(bkvalid);
    error(['An XY indicated as an alternate background is not itself',...
        ' properly defined.  Check background structure (bk) for XY(s) ',...
        num2str(xypos(~bkvalid)),'.']);
end

%Determine segmentation procedure
switch regexprep(op.local, ...
        {'(cyt.*)|(^c$)', '(nuc.*)|(^n$)'},{'cyto','nucl'})
    case 'cyto' %Cytoplasmic signals
        op.iscyto = true;
    case 'nucl' %Nuclear signals
        op.iscyto = false;
end

%Parallelization setup
if nW > 1 && matlabpool('size') == 0;
    matlabpool('local', nW);        %Start the local pool
    %Ensure availability of proper java classes from BioFormats
    if iscell(p.pth.jc); runtext = ['javaclasspath({''',p.pth.jc{1}];
        for s = 2:numel(p.pth.jc); runtext = [runtext,''',''',p.pth.jc{s}]; end%#ok<AGROW>
        runtext = [runtext,'''});'];
    else runtext = ['javaclasspath ',p.pth.jc,';'];
    end;     pctRunOnAll(runtext);
end
runid = ceil(rand*1e4);                     %Assign random run ID number


%% DATA PREPARATION  -----------------------------------------------------
%Establish Reader for ND2 file
msg = sprintf(['Loading ND2 file Reader.  ',...
    'Start time is:  ', datestr(now, 'mm/dd HH:MM PM')]);    disp(msg);
[r, ri, rp] = bfopen_custom(fname);
msg = sprintf(['Reader loaded.  ', datestr(now, 'mm/dd HH:MM PM')]);
disp(msg);
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

%Append sampling time to Global MetaData for Camera
GMD.cam.tsamp = p.tsamp;


%% PREPROCESSING PREPARATION  --------------------------------------------
%Channel ID definition (for use below)
nChan = r.getSizeC;
if nChan ~= ncname;
    warning(['Number of channels found (',num2str(nChan), ') is not ',...
        'equal to number of names provided (op.cname)']);
end

cid = cell2struct(num2cell(1:nChan), GMD.exp.Channel, 2);
%   Ensure channel reference uniformity
cn = fieldnames(cid);
%   Prepare cross talk correction if applicable
xtr = [];
if ~isempty(xtalk.D)
    xtalk = structfun(@(x) cn{ ~cellfun( 'isempty', regexpi(cn, ...
        x, 'start')) }, xtalk, 'UniformOutput', false);
    %   Get cross-talk ratios for simple removal
    xtr = iman_xtalk_correction(xtalk.I, xtalk.D, GMD, []);
end

%Save MetaData and XTalk information
save([p.sname, '_Global.mat'], 'GMD', 'xtr', 'p');
%   Short circuit if no XY positions specified (i.e. run was for meta-data)
if nxy == 0; return; end

%Determine XY storage order
if numel(ri) >= max(xypos) && ri(1).ZTCsize(1) == 1
    xySeries = true;    sz = 1;
elseif ri(1).ZTCsize(1) >= max(xypos)
    xySeries = false;   r.setSeries(0);
elseif ri(1).ZTCsize(1) == 1 && ri.ZTCsize(2) >= nChan*nxy
    xySeries = false; r.setSeries(0);
    tnxy = input(['Oversized Time dimension identified.  \nIf XY ',...
        'positions are stored in Time, provide number of XY positions:  ']);
    ri.lblr = reshape(ri.lblr, tnxy, ri.ZTCsize(2)./tnxy, nChan);
    ri.ZTCsize = size(ri.lblr);
else
    %FIXME Some sort of error I suppose
    warning(['XY storage order could not be determined. Aborting', ...
        ' and returning ND2 Reader and MetaData']);
    xtr = [];  return;
end

%Set sample display parameters  (preserves aspect ratio)
dsx = min(640,GMD.cam.PixNumX);     dsy = min(540,GMD.cam.PixNumY);
% dscx = dsx./GMD.cam.PixNumX;        dscy = dsy./GMD.cam.PixNumY;
dsc = min([dsx./GMD.cam.PixNumX, dsy./GMD.cam.PixNumY]);
dsx = GMD.cam.PixNumX*dsc;          dsy = GMD.cam.PixNumY*dsc;

%Use zero (0) as placeholder for 1 digit xypos
dhold = '0';

%Establish data output sequence (to be stored for reference)
vcorder = iman_cellmask([GMD.cam.PixNumY, GMD.cam.PixNumX, nChan], ...
    [], op, []); %#ok<NASGU>

%Get non-dynamic background intensities, if any
bki = cell(nxy,1);  hi = 1:nxy;
for h = hi(~use_altXY)  %Skips backgrounds from alternate wells
    e_inv = [];
    %     if use_altXY(h);     continue;       end
    if bk(h).fix    %IF fixed and provided
        bki{h} = bk(h).reg - bval;
    elseif ~bk(h).dyn   %IF static and region provided for image 1
        im = zeros(GMD.cam.PixNumY, GMD.cam.PixNumX, nChan);
        if xySeries; r.setSeries(xypos(h)-1);  rii = ri(p.xypos(h));
        else         sz = xypos(h);            rii = ri(1);
        end
        %Load first image for this XY
        for sc = 1:nChan; im(:,:,sc) = bfGetPlane(r,rii.lblr(sz,1,sc)); end
        %Refine image
        im = subf_refine_image(im, GMD, bval, e_inv, cid, xtalk, xtr);
        %Get background valued from first image region
        bki{h} = mean(mean( ...
            im( bk(h).reg(3):bk(h).reg(4), ...
            bk(h).reg(1):bk(h).reg(2), : ), 1 ), 2);
    end
end
%   Apply 'altxy' background mappings as needed
bki(use_altXY) = bki(altXYi(use_altXY));
%Arranged on the serial 'h' index, not the literal xypos index

%Close reader in main matlab client (will be opened in workers)
r.close();


%% Sequentially process each XY position desired
t_tot = tic;  t_perxy = [];  emsg = [];  errid = 0;  nrun = 0;
for h = [hi(~use_altXY), hi(use_altXY)] %Sequence for background containing XYs first
    t_xy = tic;         nrun = nrun + 1;
    try  %Protect other xy runs from failures
        %Ensure memory is not occupied from previous run
        clear valcube masks movieInfo tracksFinal coord
        
        %Select current Series in Reader Info, or Z slice (to index xy)
        if xySeries; rii = ri(p.xypos(h)); else sz = xypos(h); rii = ri(1); end
        %Get number of time points for this XY
            nT = rii.ZTCsize(2);
        %Determine how to slice processing
        ipw = ceil( nT./nW );               %Approx. images per Worker
        nper = ipw*ones(1,nW);  dif = nW*ipw - nT - 1; %Fix overestimate
        nper(end-dif:end) = nper(end-dif:end) - 1;  %Minus 1 excess per Worker
        
        %% Process images in parallel
        %Display status message
        t_el = toc(t_tot)/60;  t_rem = t_perxy*(nxy-h+1) - toc(t_xy)/60;
        [msg, emsg] = subf_update_timer(t_el, t_rem, msg, emsg);
        msg = sprintf(['Identifying cells for XY ', ...
            '%d (%d/%d).'], p.xypos(h), nrun, nxy);    disp(msg);
        movieInfo = cell(1,nW);  bkis = cell(nW,1); cellcount = cell(1,nW);
        parfor ps = 1:nW   %Parallel Loop
            %Load reader for each Worker
            [r2, rid2] = bfGetReader_custom(fname,false);
            if isempty(rp.memopath); r2 = loci.formats.Memoizer(r2, 100);%#ok<PFBNS>
            else r2 = loci.formats.Memoizer(r2, 100, java.io.File(rp.memopath));
            end
            r2.setId(rid2);
            
            %Select current Series in Reader, or Z slice (to index xy)
            if xySeries;  r2.setSeries(xypos(h)-1);  else r2.setSeries(0);  end %#ok<PFBNS>
            e_inv = [];  %Initialize temporary variable
            im = zeros(GMD.cam.PixNumY, GMD.cam.PixNumX, nChan); %#ok<PFBNS>
            for s = 1:nper(ps)
                    st = sum( nper(1:(ps-1)) ) + s;  %#ok<PFBNS> %Get time index
                %Load Image
                for sc = 1:nChan
                    im(:,:,sc) = double( ...    %Get current image (as double)
                        bfGetPlane( r2, rii.lblr(sz,st,sc) ) );    %#ok<PFBNS>
                end
                if dispsamp && ps==1 && s==1;   %Store orignal for display only
                    im_or{ps} = single( imresize(im, [dsy,dsx]/2) ); end
                
                %Refine Image
                im = subf_refine_image(im, GMD, bval, e_inv, cid, xtalk, xtr);
                if dispsamp && ps==1 && s==1;   %Store orignal for display only
                    im_md{ps} = single( im ); end
                
                %Get background image values
                if bk(h).dyn && ~bk(h).fix   %#ok<PFBNS>
                    if use_altXY(h);   bkit = bki{altXYi(h)}(s,:); %#ok<PFBNS>
                    else    bkit = mean(mean( ...
                            im( bk(h).reg(3):bk(h).reg(4), ...
                            bk(h).reg(1):bk(h).reg(2), : ), 1 ), 2);
                        bkis{ps}(s,:) = bkit;    %Store background used
                        if bk(h).bkonly; continue; end
                    end
                else bkit = bki{h};
                end
                %Perform Cell Idenification
                [cellcount{ps}{1,s}]= iman_cellidccount(im, op, bkit);
            end
            r2.close();
        end     %END PARALLEL LOOP
    end
    xystring = ['_xy', dhold(1:double(p.xypos(h)<10)), num2str(p.xypos(h))];
    save([p.sname, xystring, '.mat'], 'cellcount');
    %Update timing counter
    t_el = toc(t_tot)/60;   t_perxy = t_el/nrun;
    
end

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


%% SUBFUNCTIONS
%Image refinement procedure
function im = subf_refine_image(im, GMD, bval, e_inv, cid, xtalk, xtr)
%Remove baseline value from image
%   The scalar addition should be removed prior to corrections
im = iman_baseline_removal(im, bval);
%Correct objective view bias (collection efficiency)
if isempty(e_inv); [im, e_inv] = iman_objective_correction(im, GMD); %#ok<NASGU>
else    im = iman_objective_correction(im, GMD, e_inv);   end
%Correct cross-talk (optional)
if ~isempty(xtr)
    %   Perform cross-talk subtraction (here, stacked channels)
    im(:,:,cid.(xtalk.D)) = im(:,:,cid.(xtalk.D)) ...
        - xtr.(xtalk.D).(xtalk.I).*im(:,:,cid.(xtalk.I));
end
end


%Timer display update function
function [msg, emsg] = subf_update_timer(t_el, t_rem, msg, emsg)
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