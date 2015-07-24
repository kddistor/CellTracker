%CHANGEPRECISION
%   Changes numeric precsion throughout a variable.
%
%   v = changeprecision(v, p)
%   
%   changes all numeric parts of the variable v, where p specifies the new
%   precision as 'single' or 'double'.  v may be numeric, a cell array or a
%   structure.  changeprecision is called recursively to operate through
%   cells and structures.

function v = changeprecision(v, p)

if isstruct(v)           %IF v is a structure
    fn = fieldnames(v);
    %Recursively call to operate on each field, via cell extraction
    v = cell2struct( changeprecision(struct2cell(v), p), fn, 1 );
elseif iscell(v)        %IF v is a cell array
    %Recursively call to operate on each cell
    v = cellfun(@(x)changeprecision(x, p), v, 'UniformOutput', false);
elseif isnumeric(v)     %IF v is numeric
    %Perform desired precision change
    switch lower(p)
        case  'single';   v = single(v);
        case  'double';   v = double(v);
        case  'uint32';   v = uint32(v);
    end
end

end