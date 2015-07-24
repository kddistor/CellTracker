%IMAN_CELLMASK
%   Segment live cell images to generation nuclear and cytoplasmic masks,
%   and mean Cell Trace values over masked images.  Use this function after
%   a cell track has been established (e.g. via u-Track).

function [valcube, mask, gi] = iman_cellmask(imin, coord, op, bkg)

%% Operation parameters
nth = 30;           %Number of Thresholds for binary operations
dnoise = floor(0.025*op.maxD);       %Denoise size, based on expected nuc
ncgap  = floor((0.05*op.maxD)) + 2;  %Size of gap between nuc mask and cyto
ncring = floor((0.05*op.maxD)) + 1;  %Desired cyto mask thickness
ncexpand = ncgap + ncring;  %Size to expand nuclear mask for cyto

%OBSOLETE PARAMETER SELECTIONS
% ncexpand = ceil(0.10*op.maxD) + 4; 	%Size to expand nuclear mask for cyto
% ncgap = floor(0.75*ncexpand);    %Size of gap between nuclear mask and cyto
outpct = [20,80];                   %Percentiles to reject for outliers

cused = find(~isnan(bkg));  %Channel indices used
cused = reshape(cused, 1, numel(cused)); %Ensure row vector for loop use


%Scoring Values for contour filtering
maxNucArea = round(pi*op.maxD^2/4);
minNucArea = round(pi*op.minD^2/4);

%Identify input data size
if numel(imin) < 10;    sz = imin;  else    sz = size(imin);    end
if numel(sz) > 2; nchan = sz(3);  sz = sz(1:2); else nchan = 1; end	


%Ratios to calculate (numerator, denominator)
%FIXME Adjust ratio list to number of channels (for now)
% rc = regexpi(op.cname, ...
%     '(?<c>cfp|cya)?(?<y>yfp|yel)?(?<r>rfp|red)?', 'names');
r_cn = {'c','(cfp)|(cya)';'y','(yfp)|(yel)';'r','(rfp)|(red)'};
for s = 1:size(r_cn,1)
    rc.(r_cn{s,1}) = find( ~cellfun('isempty', regexpi(op.cname, r_cn{s,2})) );
end

rt = {[rc.c,rc.y], [rc.c+nchan,rc.y+nchan],...
        [rc.c,rc.r], [rc.c+nchan,rc.r+nchan]};
rt(~cellfun(@(x)length(x) == 2, rt)) = [];  nrt = numel(rt);
% rt = {[1,2], [4,5], [4,5], [1,3], [4,6], [4,6]};   nrt = numel(rt); %Obsolete
%    ^C/Y n, C/Y c, C/Y c, C/R n, C/R c  %C:CFP, Y:YFP, R:RPF c:cyto, n:nuc



%IF coordinate/track information not provided, return valcube order
%   In this case, imin may contain only the size of the image data
if isempty(coord)
    valcube = cell(2*nchan+nrt, 1);
    if ~isfield(op, 'cname') || numel(op.cname) < nchan
        %Default channel names are Ch1, Ch2, Ch3...
        op.cname = cellfun(@(x)['Ch',num2str(x)], num2cell(1:nchan),...
            'UniformOutput', false);
    end
    %Filling order is Nucs,Cytos,Ratios,Coords,Median
    lcn = {'_Nuc','_Cyto'};
    valcube(1:nchan) = cellfun(@(x)[x,lcn{1}], op.cname, ...
        'UniformOutput', false);
    valcube(nchan+(1:nchan)) = cellfun(@(x)[x,lcn{2}], op.cname, ...
        'UniformOutput', false);
    for sr = 1:nrt    %Ratios based on definitions above
        ci = mod(rt{sr}-1, nchan) + 1;  lcl = floor((rt{sr}-1)/nchan) + 1;
        valcube{2*nchan + sr} = [op.cname{ci(1)},lcn{lcl(1)},'/',...
            op.cname{ci(2)},lcn{lcl(2)}];
    end
    %Standard appendices
    valcube(end + (1:2)) = {'XCoord', 'YCoord'};    %Appended coordinates
    valcube{end + 1} = ['Note: Individual region intensities are raw.',...
        '  Ratios are background subtracted.'];
    mask = []; gi = [];
    return
end


%% Prepare for processing
%Extract segmentation channel(s) (and associated background)
imS = imin(:,:,op.chan);  imSb = bkg(op.chan);  
%Split out image channels (for indexing ease when masking)
imin = mat2cell( imin, sz(1), sz(2), ones(1,nchan) );

coord = round(coord);   %Ensure integer coordinates
%Enforce Minimum Values on coordinates (xy)
coord( (coord(:,1,2) < 1 | coord(:,1,1) < 1), :, : ) = NaN;
%Enforce Maximum Values on coordinates (xy)
coord( (coord(:,1,2) > sz(1) | coord(:,1,1) > sz(2)), :, : ) = NaN;

%Convert coordinates to linear indices, from subscripts
lind = sub2ind( sz, coord(:,:,2), coord(:,:,1) );
gi = find(lind>0);      %Store index of 'good' coordinates
%Get filtered linear index (nonzero value tracks)
lind = lind( gi );      %Note: gi is retained to fill outputs properly

%Define binary structuring elements
st1 = strel('disk', dnoise);  st2 = strel('disk', ceil(2*dnoise));
%                               ^Removes 0.5-1x dnoise size around nucleus OBSOLETE COMMENT 

%% Segment for nuclear masks
%Define intensity thresholds by quantiles
%   Only consider pixels greater than 2x background
thresholds = quantile(imS(imS > 2*imSb), linspace(0.05, 0.95, nth) );

%Invert if segmentation channel is cytoplasmic
if op.iscyto
    mxS = max(imS(:));    imS = mxS - imS;
    thresholds = mxS - thresholds;  thresholds = thresholds(end:-1:1);
    %   Thresholds oriented to always slice from low to high (re: nucleus)
    %   This is expected to limit under-masked nucleii
end

%Initialize countour matrix
ctm = zeros(sz);
%Loop to get a contour slice for each coordinate (as feasible)
tind = lind;
for s = 1:nth
    %Terminate loop if all indices identified
    if isempty(tind); break; end
    
    %Get features as thresholded
    tim = imS > thresholds(s);
    %Remove small spots (by erosion then dilation)
    tim = imdilate(tim, st1);  	%Dilate (to remove holes)
    tim = imerode(tim, st2);   	%Erode (to remove spots)
    tim = imdilate(tim, st1); 	%Redilate (to restore size, nearly)
    
    %Label connected components (disjoint features)
    tim = bwconncomp(tim);
    %Reject regions a priori if too big or not labeled by coords (tind)
    rj = cellfun(@length, tim.PixelIdxList) > maxNucArea*2;
    rj(~rj) = cellfun(@(x)~any(ismember(tind,x)), tim.PixelIdxList(~rj));
    tim.NumObjects = tim.NumObjects - sum(rj);
    %   Short circuit and continue if no valid Labels remain
    if tim.NumObjects == 0; continue; end
    tim.PixelIdxList(rj) = [];
    
    %Redo the labeling on the reduced binary image
    %Get properties for the retained features
    S = regionprops(tim, 'Area', 'Perimeter', 'EulerNumber');
    nucArea  = cat(1,S.Area);
    nucFormfactor = 4*pi*nucArea./(cat(1,S.Perimeter).^2);
    iscomplete = cat(1, S.EulerNumber) == 1;  %Check for having holes
    %Score properties (binary)
    sizeScore  = nucArea > minNucArea*.3  &  nucArea < maxNucArea;
    shapeScore = nucFormfactor > op.minF;
    totalScore = sizeScore & shapeScore & iscomplete;
    %   Rejects regions with holes (expect better segment on next slice)
    
    %Keep regions that satisfy property scores
    findex = cat(1,tim.PixelIdxList{totalScore});
    ctm( findex ) = 1;
    % Remove found contours from list of indices to search for
    tind(ismember(tind,findex)) = [];
end

%Label the nuclear region image
ctml = bwlabel(ctm);
lbl  = ctml(lind);


%% Extend to cytoplasm masks and store mean values
%Initialize valcube
ntf = size(coord,1);
valcube = nan(ntf, 1, nchan*2+nrt+2);
%Initialize mask storage structure
mask = struct('nuc', cell(ntf,1), 'cyt', cell(ntf,1));
lindex = reshape(1:numel(imS), sz);

%Generate Cytoplasm region
%Thicken nuclear region (to make a donut for cytoplasm)
%   Can simply dilate because ctml is a label matrix (different values for
%       different regions, so none will merge)
expdisk = strel('disk', ncexpand);

%Erode nuclear mask for clearance from nuc/cyt boundary
%   Use of 0.6 factor biases gap into the nucleus (qualitatively better)
ctml = imerode(ctml, strel('disk', ceil(ncgap.*0.6)));
%Expand cytoplasmic mask from final nuclear mask
cytr = imdilate(ctml, expdisk);
ctemp = ctml; ctemp(ctemp == 0) = Inf;  %Inverted dilation expands labels
ctemp = -imdilate(-ctemp, expdisk);     % in the opposite direction
ctemp(ctemp == Inf) = 0;  cytr(~(ctemp==cytr)) = 0; %Overlaps are removed

%Remove Nuclear pixels from Cytoplasm region
cytr = cytr.*~(imdilate(ctml, strel('disk', ncgap))>0);

%Calculate and fill outputs (the 'valcube' matrix)
for s = find(lbl)'
    % Get intensity values for each region.
    %Define masks
    nmask = ctml == lbl(s);
    cmask = cytr == lbl(s);
    %Store masks
    mask(gi(s)).nuc = lindex(nmask);
    mask(gi(s)).cyt = lindex(cmask);
    
    %Collect values as masked, with outliers removed (extreme 10%)
    for sc = cused
        %Nuclear mask
        vals = imin{sc}(nmask);  pctb = prctile(vals,outpct);
        valcube(gi(s), 1, sc) = mean(vals(vals > pctb(1) & vals < pctb(2)));
        %Cytoplasm mask
        vals = imin{sc}(cmask);  pctb = prctile(vals,outpct);
        valcube(gi(s), 1, nchan + sc) = ...
            mean(vals(vals > pctb(1) & vals < pctb(2)));
    end
end

%Store commonly used ratio calculations (background subtracted)
nvc = nchan*2;
for sr = 1:nrt
    bv = bkg( mod(rt{sr}-1, nchan) + 1 );               %Background values
    vct = (valcube(:,:,rt{sr}(1)) - bv(1)) ./...        %Ratio
             (valcube(:,:,rt{sr}(2)) - bv(2));          
    vct = vct.*( (valcube(:,:,rt{sr}(1)) > bv(1)) & ... %Filter out values
                 (valcube(:,:,rt{sr}(2)) > bv(2)) );    %   below zero
    valcube(:, :, nvc+sr) = vct;                        %Store
end

%Store coordinate values (to valcube slices 13 and 14)
valcube(:,:,nvc+nrt + (1:2)) = coord;

%Append a Median value cell trace
valcube(end+1, :, :) = nanmedian(valcube, 1);

%For diagnostics, show image with masking and labels
if op.diag
    figure; clf; hold on;
    nim = false(sz); cim = nim;
    glbls = unique(lbl)'; glbls(glbls==0) = [];
    for p = glbls
        nim(bwperim(ctml == p)) = true;
        cim(bwperim(cytr == p)) = true;
    end
    dimR = imS/max(imS(:)); dimG = dimR;   dimB = dimR;
    dimR(cim) = 1;   dimR(nim) = 0;
    dimG(cim|nim) = 0;
    dimB(nim) = 1;   dimB(cim) = 0;
    imshow(cat(3,dimR,dimG,dimB)); hold on;
    [yy, xx] = ind2sub(sz,lind);
    for p = find(lbl)'
        plot(xx(p),yy(p),'c.');
        text(xx(p)+5,yy(p),num2str(p),'Color','Green','FontSize',10);
    end
end

