%TRACKCOORDCLEANUP
%   Cleanup routine for extracing coordinates from tracksFinal
%   Coordinates are returned nCells x nTime x 2 (for X,Y)

function c = iman_trackcoords(tF)

%Get track coordinates from tracksFinal
nT = max( arrayfun(@(x)x.seqOfEvents(end,1), tF) ); %# of Times
nC = length( tF );                                  %# of Cells
c = nan(nC, nT, 2);    %Initialize Coordinate matrix
%FIXME compound tracks store in 1st dimention of tracksCoordAmpCG
for s = 1:nC
    %seqOfEvents gives: 1-frame index of event, 2-(1=start, 2=end),
    %   3-local track index, 
    %   4-(NaN=true start/end, value=local index of interacting track)
    truei = isnan(tF(s).seqOfEvents(:,4));      %Indices of true start/termini
    starti = tF(s).seqOfEvents(:,2) == 1;       %Indices of starts
    %TEMP: If compound tracks are reported, retain the earliest in time
    [tm, ti] = min( tF(s).seqOfEvents(starti & truei, 1) );  %#ok<ASGLU>
    li = tF(s).seqOfEvents(ti,3); %Local index (applies to other fields)
    ct1 = (tF(s).seqOfEvents(:,3) == li);   %Indices of Compound Track 1
    sf = tF(s).seqOfEvents( ct1 & truei & starti, 1);     %True start frame
    ef = tF(s).seqOfEvents( ct1 & truei & ~starti, 1);    %True end frame
    nf = ef-sf + 1; %Number of frames for Compound Track 1
    %Fill Coordinates from tracksFinal (Time lengths may vary)
    c(s, sf:ef, 1) = tF(s).tracksCoordAmpCG(li,1:8:8*nf); %x-Coord
    c(s, sf:ef, 2) = tF(s).tracksCoordAmpCG(li,2:8:8*nf); %y-Coord
end
%NaNify zero values coordinates
c(c==0) = NaN;

%Interpolate gaps in coordinate vectors
%   Fill bounded NaNs with interpolated values
index = 1:size(c,2);
for s = 1:nC
    %Get failed coordinates
    nanD = any( isnan(c(s,:,:)), 3 );
    if ~any(nanD); continue; end    %Skip if no failed coordinates
    %Skip any track with < 2 valid coordinates
    if nnz(~nanD) < 2; c(s, ~nanD, :) = NaN;  continue;   end
    %Interpolate between valid coordinates
    c(s, nanD, :) = interp1(index(~nanD), shiftdim(c(s, ~nanD, :),1), ...
        index(nanD), 'linear', NaN); 
end


end