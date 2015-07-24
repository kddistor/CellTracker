%IMAN_CELLID
%   Image Analysis Cell Identification procedure.

function [movieInfo, cellmasks] = iman_cellid2(im, op, bkg)

nth = 20;           %Number of Thresholds for binary operations

%% Preliminary calculations
%Calculations used for cell sizing
maxNucArea = round(pi*op.maxD^2/4);
minNucArea = round(pi*op.minD^2/4); %0.3 from below removed
flts = floor(op.maxD/10);           %Filter size, based on expected nuc

%Background removal, if supplied
if exist('bkg', 'var') && ~isempty(bkg)
    im = bsxfun( @minus, im, reshape(bkg, 1, 1, numel(bkg)) );
    ebkg = bkg(op.chan);
else
    ebkg = 0;
end

%Define a Gaussian filter (smoothes noise)
gaussianFilter = fspecial('gaussian', flts.*[1, 1], 10);
e = double( imfilter(im(:,:,op.chan), gaussianFilter, 'replicate') );
e_for = e > ebkg;   %Foreground area mask

%Define binary structuring elements
st1 = strel('disk', ceil(flts/4));  st2 = strel('disk', 2*ceil(flts/4));


%% Image filtering (IF necessary, i.e. using edge filters for cyto)
if op.iscyto
    %Create a predefined 2d filter. "Sobel edge-emphasizing filter"
    hy = fspecial('sobel');   hx = hy';
    %Filter image for edges (gradient, via the Sobel filter)
    e = sqrt(imfilter(e, hx, 'replicate').^2 ... %Across x axis.
             + imfilter(e, hy, 'replicate').^2); %Across y axis.
    e = max(e(:)) - e;                  %Invert to make gradients 'low' 
end

%% Image thresholding search
%Threshold at linearly spaced quantiles
thresholds = quantile(e(e_for), linspace(0.05, 0.95, nth) ); %Thresholds
clear e_for;

%FOR each threshold, evalute appearance of nucleus-like shapes
enuc = false(size(e));
for j=1:length(thresholds)
    %Threshold the gradient image and erode for noise
    et = e>thresholds(j);           %Threshold
    et = imdilate(et, st1);       	%Dilate (to remove holes)
    et = imerode(et, st2);         	%Erode (to remove spots)
    et = imdilate(et, st1);         %Redilate (to restore size)
    
    %Segment and label the binary image
    etl = bwlabel(et);
    %Get properties of each labelled region (Area, Perimeter)
    S = regionprops(etl, 'Area', 'Perimeter'); %#ok<MRPBW>
    nucArea  = cat(1,S.Area);                      %Vector of Areas
    nucPerim = cat(1,S.Perimeter);                 %Vector of Perimeters
    nucFormfactor = 4*pi*nucArea./(nucPerim.^2);   %Area/Perim ratio
    
    %Calculate binary 'scores' from Area and form factor (Perimeter)
    %   Must be between min and max Areas
    sizeScore  = nucArea > minNucArea  &  nucArea < maxNucArea;        
    %   Must be greater than min (ratio of Area to Perimiter like a circle)
    shapeScore = nucFormfactor > op.minF;              
    
    %Get combined score (binary indicator of satisfaction)
    totalScore = sizeScore & shapeScore;
    scorenz = find(totalScore);
%     etl=imerode(etl, strel('disk',4));
    %Filter by score and store estimated cell nuclei
    enuc( ismember(etl, scorenz) ) = true;
    
end

% enuc = imerode(enuc, strel('disk',3));


%% Store estimate nuclei data
%Re-segment stored nuclei spot (re-orders spots)
C = regionprops(enuc, 'Centroid', 'Area', 'PixelIdxList');
nCtr = cat(1, C.Centroid);
nAre = cat(1, C.Area);
im_op = im(:,:,op.chan);
nAmp = arrayfun(@(x)sum(im_op(x.PixelIdxList)), C);
cellmasks = bwconncomp(enuc);
%Store coordinates IF any spots were found
if ~isempty(nCtr)
    movieInfo.xCoord = nCtr(:,1);
    movieInfo.yCoord = nCtr(:,2);
    movieInfo.amp = nAmp(:,1);
    movieInfo.are = nAre(:,1);
else  %Store empty output if no spots found
    movieInfo = struct('xCoord',[],'yCoord',[],'amp',[], 'are',[]);
end

