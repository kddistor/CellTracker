%BFOPEN_CUSTOM
%   Opens ND2 files via the BioFormats toolbox for Matlab, with option to
%   avoid loading all images at once (to reduce memory load).
%
%   [r, out] = bfopen_custom(in, varargin)
%
%   returns the Java reader object, r, which links to the ND2 file, and an
%   output structure containing labeling data, metadata, and images (if
%   requested).
%
%   out structure contains fields
%       im:    all image data from the ND2 file, only if bfopen_custom is
%               called with the additional input pair ('all', true)
%       ZTCsize:    data size in order Z (often XY), Time, Channel
%       lbl:   label information for XY, Time, and Channel, organized as a
%               structure, with fields Z, T, C (Z  may represent XY index)
%       lblr:  reverse label information, a 3-D matrix which gives the
%               Plane index (for accessing images via the Reader, r) by the
%               Z, T, C coordinates.  lblr(Z,T,C) gives the Plane index.
%       time:   Timestamp values
%
%Note:  Ensure the file is local.  Data must be pulled at some point and if
%   remote, the network limitations will take effect when parsing the file
%   (r.setId in bfGetReader)


function [r, out, p] = bfopen_custom(in, varargin)

%Default options
p.all = false;  %Extract and return all images
% p.memopath = '\\cloud.mcb.ucdavis.edu\albeck_lab\imageData\bf_memo';
p.memopath = [];

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end %#ok<WNTAG>
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%Initialize reader object (Time-consuming step)
[r, rid] = bfGetReader_custom(in,false);
%Decorate with 'Memoizer' - to cache the reader for large files, such that
%   initialization will be faster in the future
if isempty(p.memopath);     r = loci.formats.Memoizer(r, 100);
else    r = loci.formats.Memoizer(r, 100, java.io.File(p.memopath));
end
%Initialize reader to the desired file, caching as needed
r.setId(rid);

%Get series count
nser = r.getSeriesCount();
%   Series must be done in loop, using r.setSeries()

out = struct('im', [], 'ZTCsize', [], 'lbl', [], 'lblr', [], 'time', []);
for s = 1:nser  %FOR Each Series
    %% Get labeling (ZCT coords) and images if requested
    %Set Reader to the current Series
    r.setSeries(s-1);
    
    %Get image count
    nim = r.getImageCount();
    %Collect images if getting all at once
    if p.all
        out(s).im = arrayfun(@(x)bfGetPlane(r, x), (1:nim), ...
            'UniformOutput', false);
    end
    
    %Collect labeling information (Z, Time, Channel indices)
    ltemp = arrayfun(@(x)label_info(r,x), (1:nim), 'UniformOutput', false);
    lbl = [ltemp{:}];  clear ltemp;
    %Prepare reverse index for Plane number by Z, Time, Channel
    lblr = zeros(lbl(end).Z, lbl(end).T, lbl(end).C);
    lblr(sub2ind([lbl(end).Z, lbl(end).T, lbl(end).C], ...
                [lbl.Z], [lbl.T], [lbl.C])) = 1:numel(lbl);
    out(s).lbl = lbl;  out(s).lblr = lblr;
    out(s).ZTCsize = [r.getSizeZ,r.getSizeT,r.getSizeC];
    
    
    %% Get time values from series metadata
    %Get series metadata
    smd = r.getSeriesMetadata();
    %Collect MetaData names
    n = arrayfun(@(x)char(x), smd.keySet.toArray, 'UniformOutput', false);
    %Collect MetaData values
    v = arrayfun(@(x)x, smd.values.toArray, 'UniformOutput', false);
    
    %Sort for timestamps and times (in case of other entries)
    %   Get Index values for 'timestamps'
    ti = cellfun(@(x)regexpi(x, '^timestamp\s+#?(?<tval>\d+)',...
            'tokens'), n, 'UniformOutput', false);
    %   Get logical to extract only timestamp data
    tivalid = ~cellfun('isempty', ti);
    %   Get all timestamps and convert to double
    ti = cellfun(@(x)str2double(x{1}{1}), ti(tivalid)); 
    %   Check for zero valued time point and correct
    if min(ti) == 0; ti = ti + 1; end
    
    %Map time values to time indices
    tv = cat(1,v{tivalid});  out(s).time(ti,1) = tv;
        
end  %FOR Each Series

end


%Label Info function
function [linfo] = label_info(rin, x)
z = num2cell(rin.getZCTCoords(x-1) + 1);    %Get Position, Color, Time
%   Adjusted for use of 0 as first coordinate
[linfo.Z, linfo.C, linfo.T] = deal(z{:});   %Deal to Structure
end


