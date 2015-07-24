%IQ_GETSPECTRALPAR
%   Gathers parameters for image quantification from stored data files.
%
%   iqp = iq_getspectralpar()
%       returns the structure iqp with all Fluorophore and Filter spectra
%       and parameters

function iqp = iq_getspectralpar()
%Wavelength range
wl_lo = 300;  wl_hi = 800;      %nm
iqp.WaveLength = wl_lo:wl_hi;   %Expect always 1 nm step

%Define data files to use
data_properties = ['\\mcb.ucdavis.edu\Shared\Data\mcb_labs\albeck\',...
    'Notebooks\MichaelPargett\WorkingShare\Imaging Quantification\',...
    'Filter and fluorophore data\',...
    'FluorophoreProps_Compiled_2014-04-08.xlsx'];
data_spectra = ['\\mcb.ucdavis.edu\Shared\Data\mcb_labs\albeck\',...
    'Notebooks\MichaelPargett\WorkingShare\Imaging Quantification\',...
    'Filter and fluorophore data\',...
    'FluorophoreSpectra_Compiled_2014-04-08.xlsx'];
data_filters = ['\\mcb.ucdavis.edu\Shared\Data\mcb_labs\albeck\',...
    'Notebooks\MichaelPargett\WorkingShare\Imaging Quantification\',...
    'Filter and fluorophore data\',...
    'FilterSpectra_Compiled_2014-04-08.xlsx'];

%Get Filter/Fluorophore naming conventions
[fpn, ftn] = iman_naming();

%Get data
[dprop.n, dprop.t, dprop.r] = xlsread(data_properties);
[dspec.n, dspec.t, dspec.r] = xlsread(data_spectra);
[dfilt.n, dfilt.t, dfilt.r] = xlsread(data_filters);

%% PARSE FOR HEADERS
%Main header index for data in properties file
dprop.mainhead = find( ~cellfun('isempty', ...
    regexpi(dprop.t(:,1), '^(fluoro)(phore|chrome)?')) );
%Standard Deviation header index for properties file
dprop.sdhead = find( ~cellfun('isempty', ...
    regexpi(dprop.t(:,1), '^(sd|(st|standard)\s?(dev|deviation))')) );

%Main header index for data in spectra file
dspec.mainhead = find( ~cellfun('isempty', ...
    regexpi(dspec.t(:,1), '^(fluoro)(phore|chrome)?')) );

%Main header index for data in spectra file
dfilt.mainhead = find( ~cellfun('isempty', ...
    regexpi(dfilt.t(:,1), '^app(lication)?')) );


%% COLLECT FLUOROPHORE PROPERTIES
%Get names used in the properties data file
names_in = lower( dprop.t(dprop.mainhead,:) );
f_ind = cellfun(...                     %Get sheet index for each name
    @(x)~cellfun('isempty', regexpi(names_in, x)),...
        fpn(:,2), 'UniformOutput', false);
    f_use = cellfun(@(x)any(x), f_ind);     %Remove unmatched names
fmatch = cell2struct(f_ind(f_use), fpn(f_use,1), 1); %Structure with index

%Fill values for each fluorophore
ind_QY = dprop.mainhead + find(~cellfun('isempty', ...
    regexpi(dprop.t(dprop.mainhead+1:dprop.sdhead,1), '(qy|quantum\syield)') ));
ind_mec = dprop.mainhead + find(~cellfun('isempty', ...
    regexpi(dprop.t(dprop.mainhead+1:dprop.sdhead,1), ...
        '(mec|m(olar)?.+ext(inction)?)') ));
for s = find(f_use)'  %Fill fluorophore properties from file
    iqp.(fpn{s,1}).QY = dprop.r{ind_QY, fmatch.(fpn{s,1})}; 
    iqp.(fpn{s,1}).mec = dprop.r{ind_mec, fmatch.(fpn{s,1})};
end


%% FILL FLUOROPHORE SPECTRA
%   For each name, get indices from the sheet.  Then get ex/em.  Then take
%   wavelengths from 300:800.  Convert NaN to 0.  Store as column vectors.
%   
names_in = lower( dspec.t(dspec.mainhead,:) );  %Names from sheet
f_ind = cellfun(...                     %Get sheet index for each name
    @(x)~cellfun('isempty', regexpi(names_in, x)),...
        fpn(:,2), 'UniformOutput', false);
f_use = cellfun(@(x)any(x), f_ind);     %Remove unmatched names
fmatch = cell2struct(f_ind(f_use), fpn(f_use,1), 1); %Structure with index

typen = {'Ex', 'ex(itation)?'; 'Em', 'em(ission)?'};
ind_type = ~cellfun('isempty',regexpi(dspec.t(:,1), '^type'));
%Get indices for each type of spectra
for s = 1:size(typen,1)
    ftype.(typen{s,1}) = ~cellfun('isempty',regexpi(dspec.t(ind_type,:), typen{s,2}));
end

%Get indices for the wavelengths (300 to 800 nm, incremented by 1)
ind_wl = dspec.mainhead + find(~cellfun(...
    'isempty',regexpi(dspec.t(dspec.mainhead+1:end,1), 'wave(length)?') ));
ind_wl = ind_wl + ( find([dspec.r{ind_wl+1:end,1}] == wl_lo) : ...
                    find([dspec.r{ind_wl+1:end,1}] == wl_hi) );
%Finally, take spectral data for each, stored by uniform name
for s = find(f_use)'
    iqp.(fpn{s,1}).ab = [dspec.r{ind_wl, fmatch.(fpn{s,1}) & ftype.Ex}];
    iqp.(fpn{s,1}).em = [dspec.r{ind_wl, fmatch.(fpn{s,1}) & ftype.Em}];
end


%% FILL FILTER SPECTRA
%   For each name, get indices from the sheet.  Then get ex/em.  Then take
%   wavelengths from 300:800.  Convert NaN to 0.  Store as column vectors.
%   
names_in = lower( dfilt.t(dfilt.mainhead,:) );  %Names from sheet
f_ind = cellfun(...                     %Get sheet index for each name
    @(x)~cellfun('isempty', regexpi(names_in, x)),...
        ftn(:,2), 'UniformOutput', false);
f_use = cellfun(@(x)any(x), f_ind);     %Remove unmatched names
fmatch = cell2struct(f_ind(f_use), ftn(f_use,1), 1); %Structure with index

typen = {'Ex', 'ex(itation)?'; 'Em', 'em(ission)?'; 'Mr', '(dichroic)?\s?mirror'};
ind_type = ~cellfun('isempty',regexpi(dfilt.t(:,1), '^type'));
%Get indices for each type of spectra
for s = 1:size(typen,1)
    ftype.(typen{s,1}) = ~cellfun('isempty',regexpi(dfilt.t(ind_type,:), typen{s,2}));
end

%Get indices for the wavelengths (300 to 800 nm, incremented by 1)
ind_wl = dfilt.mainhead + find(~cellfun(...
    'isempty',regexpi(dfilt.t(dfilt.mainhead+1:end,1), 'wave(length)?') ));
ind_wl = ind_wl + ( find([dfilt.r{ind_wl+1:end,1}] == wl_lo) : ...
                    find([dfilt.r{ind_wl+1:end,1}] == wl_hi) );
%Finally, take spectral data for each, stored by uniform name
for s = find(f_use)'
    %Mirror refelcts/transmits light dependent on wavelength.  
    %Take product with Mirror for Em, and with 1-Mirror for Ex.
    iqp.(ftn{s,1}).ex = [dfilt.r{ind_wl, fmatch.(ftn{s,1}) & ftype.Ex}]...
        .*(1-[dfilt.r{ind_wl, fmatch.(ftn{s,1}) & ftype.Mr}]);
    iqp.(ftn{s,1}).em = [dfilt.r{ind_wl, fmatch.(ftn{s,1}) & ftype.Em}]...
        .*[dfilt.r{ind_wl, fmatch.(ftn{s,1}) & ftype.Mr}];
end

%Remove all NaNs, replacing with zeros
iqp = zeronan(iqp);

iqp.FPhoreNames = fpn;     iqp.FilterNames = ftn;

%Copy identical spectra, as needed
iqp.mCer3.ab = iqp.mCer.ab;     iqp.mCer3.em = iqp.mCer.em;

%IF no outputs requested, save the output in a MAT file
if nargout == 0;	save('iq_spectral', 'iqp');    end

end


%Sub-Function to set NaN to zero
function out = zeronan(in)
switch class(in)
    case 'struct'; out = structfun(@zeronan, in, 'UniformOutput', false);
    case 'cell';   out = cellfun(@zeronan, in, 'UniformOutput', false);
    case 'double'; out = in; out(isnan(in)) = 0;
    otherwise;     out = in; 
end
end