%IMAN_NAMING
%   Naming conventions for Filters and Fluorophores
%   
%   [fpn, ftn] = iman_naming()
%       To get cell arrays with current naming conventions and regular
%       expression maps.

function [fpn, ftn] = iman_naming()

%Define Fluorophore names, to support consistency
fpn = { 'DAPI',         'dapi'; ...
        'SECFP',        's?e\s?cfp'; ...
        'mCer',         'm?cer(ulean)?\D*$'; ...
        'mCer3',        'm?cer(ulean)?\s*3\D*$'; ...
        'mTurq2',       'm?turq(uoise)?\s*2\D*$'; ...
        'TSap',         't?.?sap(phire)?'; ...
        'YPet',         'y\s?pet'; ...
        'EGFP',         'e\s?gfp'; ...
        'Emerald',      'emer(ald)?'; ...
        'Clover',       'clov(er)?'; ...
        'mVenus',       'm?ven(us)?'; ...
        'mCherry',      'm?che?(rry)?'; ...
        'Alex647',      'alexa?(fluor)?\s*647\D*$'; ...
        'Alex680',      'alexa?(fluor)?\s*680\D*$'          };

%Define Filter names, to support consistency
ftn = { 'Filter_DAPI',         'dapi'; ...
        'Filter_CFP',          '(49001)|cfp'; ...
        'Filter_Sapphire',     'sap(phire)?'; ...
        'Filter_GFP',          'gfp'; ...
        'Filter_YFP',          '(49003)|yfp'; ...
        'Filter_Cherry',       '(49005)|dsred|dia|cy3|rfp|che?(rry)?(\s|_|-){0,3}$'; ...
        'Filter_Cherry2',      'che?(rry)?(\s|_|-){0,3}2|tritc'; ... TEMPORARY 'tritc'
        'Filter_Cy5',          'cy\s?5'           };

end