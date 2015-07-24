%ND2_META
%   Extract global meta data from ND2 file
%

%NOT YET COMPLETE.  Update help header with info for expected field values.


function md = nd2_meta(r,ri,bkin)
%Metadata structure design
%      FieldName         Text Search           Search Location  Units 
%Objective Metadata
mdi.obj.Desc =           {'.+plan.+apo.+',           'value',    ''};
mdi.obj.NA =             {'numerical\s+aperture',    'key',      'none'};
mdi.obj.Mag =            {'magnification',           'key',      'none'};
mdi.obj.WkDist =         {'working\s+dist',          'key',      'mm'};
mdi.obj.RefIndex =       {'refract\w+\s+\index',     'key',      'none'};

%Camera Metadata
mdi.cam.Desc =           {'camera\s+name',           'key',      ''};
mdi.cam.PixSizeX =       {'calibration',             'key',      'um'};
mdi.cam.PixSizeY =       {'calibration',             'key',      'um'};
mdi.cam.PixNumX =        {'uiwidth$',                'key',      'none'};
mdi.cam.PixNumY =        {'uiheight$',               'key',      'none'};
mdi.cam.BinSizeX =       {'dbinningx',               'key',      'none'};
mdi.cam.BinSizeY =       {'dbinningy',               'key',      'none'};

%Experiment Metadata may be multi-valued, per Channel
nc = r.getSizeC;
for s = 1:nc
mdi.exp.Channel(s,:) =   {['name\s+#',num2str(s)],           'key',      ''};
mdi.exp.Exposure(s,:) =  {['exposure\s+#',num2str(s)],       'key',      'ms'};
mdi.exp.ExVolt(s,:) =    {['sola.+voltage\s+#',num2str(s)],  'key',      ''};
mdi.exp.Filter(s,:) =    {['filter.+turret.+\s+#',num2str(s)],  'key',   ''};
mdi.exp.FPhore(s,:) =    {['fluorophore',num2str(s)],            'key',  ''};
end

    
%% Get Global Metadata
gmd = r.getGlobalMetadata();
%   Collect MetaData names
n = arrayfun(@(x)char(x), gmd.keySet.toArray, 'UniformOutput', false);
%   Collect MetaData values
v = arrayfun(@(x)x, gmd.values.toArray, 'UniformOutput', false);

%Reduce Metadata size
gind = cellfun( @(x)isempty(regexpi(x, '^(X|Y|Z|time)', 'start')), n );
n = n(gind); v = v(gind);

%Extract Metadata to the defined format
md = extract_mdata(n, v, mdi);

%Display Metadata from ND2 file
display('Showing Metadata extracted from the ND2 file:');
for s = fieldnames(md)';  display(s{1});  display(md.(s{1}));   end

%Check defined Metadata for empty fields and request values
if ~exist('bkin','var'); bkin = []; else display('Using backup.'); end
md = check_mdata_complete(md, bkin);

%Metadata clean up (ensure naming uniformity, unit nomenclature)
md = md_cleanup(md, mdi);

%Double-check Metadata completeness (Clean up may clear fields)
md = check_mdata_complete(md, bkin);

end



%% Metadata extraction subfunction (recursive on structure)
function out = extract_mdata(keys, vals, wants)

if isstruct(wants)
    %IF structure, recursively fill each field
    fnames = fieldnames(wants);
    for s = 1:length(fnames)
        out.(fnames{s}) = extract_mdata(keys, vals, wants.(fnames{s}));
    end
else
    %If not a structure, prepare the output and extract
    out = cell(1,size(wants,1));
    for s = 1:size(wants,1)
        switch lower(wants{s,2})
        case 'key'  %IF searching on the keys (field names)
            try %#ok<TRYNC>
            out(1,s) = vals( cellfun( ...
                @(x)~isempty(regexpi(x, wants{s,1}, 'start')), keys ) );
            end
        case 'value'  %IF searching on the values
            try %#ok<TRYNC>
            out(1,s) = vals( cellfun( ...
                @(x)~isempty(regexpi(x, wants{s,1}, 'start')), ...
                cellfun(@(x)num2str(x),vals,'UniformOutput',false) ) );
            end
        end
    end
        %Reduce cell array if feasible, converting to doubles
        %   Try converting to doubles   
        temp = cellfun(@(x)str2double(x), out);
        %   Only use if all cells converted properly
        if ~any(isnan(temp))
            out = temp;
        %Reduce cell if only 1 value or all values are numeric
        elseif ~any(cellfun('isempty',out(1,s))) && ...
            ( length(out) == 1 || all(cellfun(@isnumeric,out)) )
            out = cell2mat(out);
        end
        %Enforce positive numeric values
        if isnumeric(out);  out(out < 0) = [];  end
        
end
end



%% Metadata completeness checking
function md = check_mdata_complete(md, bkin, prnt)
%Initialize
if ~exist('prnt','var'); prnt = '';
elseif isstruct(bkin) && isfield(bkin,prnt); bkin = bkin.(prnt); end 
isfirstwarn = true;

%Checking routine
fn = fieldnames(md);
for s = 1:length(fn)
    if isstruct(md.(fn{s}))
        %IF structure, recursively check each field
        if isempty(prnt); spc = ''; else spc = '.'; end
        md.(fn{s}) = check_mdata_complete(md.(fn{s}), bkin, [prnt,spc,fn{s}]);
    else
        %If not a structure, check / request values / return
        if iscell(md.(fn{s}))
            unpop = any(cellfun('isempty',md.(fn{s})));
        else
            unpop = isempty(md.(fn{s}));
        end
        
        if unpop
            %First check for data in the backup input
            try %#ok<TRYNC>
                md.(fn{s}) = bkin.(fn{s});
                continue;
                %If found, go to next loop
            end
            %If First unpopulated field, print warning info
            if isfirstwarn
                warning(['Unpopulated Metadata found. Enter values as ',...
                    'requested. Refer to "help nd2_meta" for field ',...
                    'descriptions.']);
                isfirstwarn = false;
            end
            %   Get current value of bad field
            %Request input to fill field
            md.(fn{s}) = input(['\nEnter value(s) for: ', prnt,'.',fn{s}, ...
                '\n    This field ',...
                'requires ', num2str(length(md.(fn{s}))), ' value(s) ',...
                '\n    New value(s): ']);
        end
    end
end

end


%% Metadata cleanup
function md = md_cleanup(md, mdi)

%Uniform naming conventions
if isfield(md,'exp')
    [fpn, ftn] = iman_naming();  %Get naming convention
    %   Fluorophore naming
    for s = 1:numel(md.exp.FPhore)
        fmatch = ~cellfun('isempty', ...
            regexpi(md.exp.FPhore{s},fpn(:,2),'start'));
        %Ensure a match is found or eliminate field
        if ~any(fmatch)
            warning(['Fluorophore "',md.exp.FPhore{s},'" not found.',...
                '  Clearing field for backup check.']);
            md.exp.FPhore{s} = [];
        else
            md.exp.FPhore{s} = fpn{fmatch,1};
        end
    end
    %   Filter naming
    for s = 1:numel(md.exp.Filter)
        fmatch = ~cellfun('isempty', ...
            regexpi(md.exp.Filter{s},ftn(:,2),'start'));
        %Ensure a match is found or eliminate field
        if ~any(fmatch)
            warning(['Filter "',md.exp.Filter{s},'" not found.',...
                '  Clearing field for backup check.']);
            md.exp.Filter{s} = [];
        else
            md.exp.Filter{s} = ftn{fmatch,1};
        end        
    end
    %   Channel naming (for valid fieldnames)
    md.exp.Channel = cellfun(@(x)regexprep(x, '\W|^\d', '_'), ...
        md.exp.Channel, 'UniformOutput', false);
end


%Numeric value and unit extraction
fn = fieldnames(md);
for s = 1:numel(fn)
   if isstruct(md.(fn{s}))  %Recursive through structure (again)
       md.(fn{s}) = md_cleanup( md.(fn{s}), mdi.(fn{s}) );
   else
       if ~isempty(mdi.(fn{s}){1,3}) && ~isnumeric(md.(fn{s}))
       %IF a unit is requested and input is a string or cell array
       switch lower(mdi.(fn{s}){1,3})
           case 'none'  %In Case of unitless result (want numeric)
           val = regexpi(md.(fn{s}),'\d+\.?\d+','match');
           md.(fn{s}) = str2double(val{1});
           otherwise    %In case units requested (ID units provided)
           val = regexpi(cellfun(@num2str, md.(fn{s}), 'UniformOutput', false),...
               '(?<num>\d+\.?\d+)\s*(?<unt>\w+)?','names');
           md.(fn{s}) = unit_chk([val{:}], mdi.(fn{s}){1,3});
       end
       end       
   end
end
    
end


function val = unit_chk(in, tun)
pfix = struct('p',1,'n',2,'u',3,'m',4,'k',6,'M',7);

%Convert to double
val = arrayfun(@(x)str2double(x.num), in);

%Check units included
mtch = regexp({in.unt}, ...
    '(?<factor>p|n|u|m|c|d|k|M?)\s?(?<unit>s|sec|m|%?)', 'names');
mtcht = regexp(tun, ...
    '(?<factor>p|n|u|m|c|d|k|M?)\s?(?<unit>s|sec|m?)', 'names');

mtch = [mtch{:}]; %Removes any empty results (from a unitless value)

%Modify value to match target units, as needed
if ~isempty(mtch(1).unit) && ~isempty(mtch(1).factor) %Uses first good unit
    val = val.*10^(3*( pfix.(mtch(1).factor) - pfix.(mtcht.factor) ));
end

end

%%

%Notes for extraction from OMEXMLMetadataStore
%   http://downloads.openmicroscopy.org/bio-formats/5.0.3/api/ome/xml/meta/OMEXMLMetadataImpl.html
%   getPixelsPhysicalSizeX(int imageIndex) 
%   getPixelsPhysicalSizeY(int imageIndex) 
%   getPixelsPhysicalSizeZ(int imageIndex) 
%   getPixelsSizeC(int imageIndex) 
%   getPixelsSizeT(int imageIndex) 
%   getPixelsSizeX(int imageIndex) 
%   getPixelsSizeY(int imageIndex) 
%   getPixelsSizeZ(int imageIndex) 
%   getPlaneDeltaT(int imageIndex, int planeIndex) 
%   getPlanePositionX(int imageIndex, int planeIndex)
%   getPlanePositionY(int imageIndex, int planeIndex)
%   getPlanePositionZ(int imageIndex, int planeIndex)
%   getObjectiveLensNA(int instrumentIndex, int objectiveIndex) 
%   getObjectiveSettingsRefractiveIndex(int imageIndex) 
%   getObjectiveWorkingDistance(int instrumentIndex, int objectiveIndex) 
%   getObjectiveNominalMagnification(int instrumentIndex, int objectiveIndex) 


% %Display image
% figure();
% for s = [1:12]
% subplot(4,3,s);
% if s > 9; s = s + 63; end
% imagesc(z.data{s});
% title(['Plane: ', num2str(s),'; XYTind: ', num2str(z.Z{s}),...
%     '; Channel: ', num2str(z.C{s})]);
% end

