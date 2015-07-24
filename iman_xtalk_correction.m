%IMAN_XTALK_CORRECTION
%
%   Image, MetaData, and Parameters (ind/dep signals) must be provided.
%	Image must have DC offset removed (the 100 from the camera, e.g.)

function xtr = iman_xtalk_correction(Ind, Dep, MD, SP)

%% Input management
%Ensure input consistency (as cell) for Channel Flags
if ~iscell(Ind); Ind = {Ind}; end; if ~iscell(Dep); Dep = {Dep}; end 

%Ensure spectral parameters are available
if ~exist('SP','var') || isempty(SP)
    %Run routine to extract spectral parameters
    SP = iq_getspectralpar;
end

%Get excitation power estimates
exPower = iq_excitation_powermap(MD.exp.Filter, MD.exp.ExVolt);

%Get camera QE (quantum efficiency) curve
QEcam = iq_getcameraqe(MD.cam.Desc, SP.WaveLength);

nd = numel(Dep);  ni = numel(Ind);

%Collect indices per Channel (to access Filter/FPhore spectra)
%   For Dependent Channels
for s = 1:nd
    iD(s) = find( ~cellfun('isempty', ...
                regexpi(MD.exp.Channel, Dep{s}, 'start')) ); %#ok<AGROW>
end
%   For Independent Channels
for s = 1:ni
    iI(s) = find( ~cellfun('isempty', ...
                regexpi(MD.exp.Channel, Ind{s}, 'start')) ); %#ok<AGROW>
end


%% X-Talk Ratio estimation
%Fundamental procedure:
%   Dependent Image = Dependent Intensity + sum(Independent Contributions)
%       Independent Contribution = Independent Image .* X-Talk_Ratio
%       Ratio is photons to Dep. Image ./ photons to Ind. Image;
%
%   Independent channel intensity = 
%       C_ind .* sum(ExF_ind./sum(ExF_ind).*ExS_ind./sum(ExS_ind)) .* dL .* MEC_ind
%             .* sum(EmF_ind.*QE_cam.*EmS_ind./sum(EmS_ind)) .* tExp_ind
%             .* QY_ind .* eta .* Gain .* exPower_ind
%
%   Dependent channel contribution = 
%       C_ind .* sum(ExF_dep./sum(ExF_dep).*ExS_ind./sum(ExS_ind)) .* dL .* MEC_ind
%             .* sum(EmF_dep.*QE_cam.*EmS_ind./sum(EmS_ind)) .* tExp_dep
%             .* QY_ind .* eta .* Gain .* exPower_dep
%
%   Non-spectral terms cancel in the X-Talk Ratio =
%       sum(ExF_dep.*ExS_ind) .* sum(EmF_dep.*QE_cam.*EmS_ind) 
%       .* tExp_dep.*exPower_dep  ./  
%     ( sum(ExF_ind.*ExS_ind) .* sum(EmF_ind.*QE_cam.*EmS_ind) 
%       .* tExp_ind.*exPower_ind )


%Calculate X-Talk Ratios (For each Dependent/Independent pair)
for sd = 1:nd
    for si = 1:ni
    %X-Talk Ratio evaluated as the product of three terms (each a ratio):
    xtr.(MD.exp.Channel{iD(sd)}).(MD.exp.Channel{iI(si)}) = ...
      ...%Total Photon Delivery (per Wavelength bin)
    (MD.exp.Exposure(iD(sd)) .* ...                 %Dependent
        exPower.(MD.exp.Filter{iD(sd)}) ./ ...      
        sum(SP.(MD.exp.Filter{iD(sd)}).ex) )        ./      ...  
    (MD.exp.Exposure(iI(si)) .* ...                 %Independent
        exPower.(MD.exp.Filter{iI(si)}) ./ ...      
        sum(SP.(MD.exp.Filter{iI(si)}).ex) )    ... 
    .*...%Excitation Rate
    sum( SP.(MD.exp.Filter{iD(sd)}).ex .* ...       %Dependent
         SP.(MD.exp.FPhore{iI(si)}).ab )            ./      ...    
    sum( SP.(MD.exp.Filter{iI(si)}).ex .* ...       %Independent
         SP.(MD.exp.FPhore{iI(si)}).ab + eps ) ...
    .*...%Emission/Collection Rate
    sum( SP.(MD.exp.Filter{iD(sd)}).em .* ...       %Dependent
         SP.(MD.exp.FPhore{iI(si)}).em .* QEcam )   ./      ...
    sum( SP.(MD.exp.Filter{iI(si)}).em .* ...       %Independent
         SP.(MD.exp.FPhore{iI(si)}).em .* QEcam + eps );%eps to prevent NaN
    end
end



%% Notes on full quantification

% %Need estimate of illuminated area
% Gain = 1;  Ceff = 0.17;
% 
% %Non-linear version
% wl_ex = sum(SP.WaveLength .* SP.(MD.exp.Filter{iI(si)}).ex...
%     ./ sum( SP.(MD.exp.Filter{iI(si)}).ex ));    %Mean excitation wavelength
% 
% l_path = 10^-4*...
%     ( (10^-3) .* wl_ex .* MD.obj.RefIndex ./ (MD.obj.NA.^2))...
%     + ( MD.obj.RefIndex .* MD.cam.PixSizeX ./ (MD.obj.Mag .* MD.obj.NA) );
% 
% %Quantitative estimate
% cI = log10( 1 - imI./(Gain.*Ceff.*SP.(MD.exp.FPhore{iI(si)}).QY.*...
%     MD.exp.Exposure(iI(si)).*exPower.(MD.exp.Filter{iI(si)}) .* ...
%     QEcam.* ...
%     SP.(MD.exp.Filter{iI(si)}).em .* ...
%     SP.(MD.exp.FPhore{iI(si)}).em .* ...
%     SP.(MD.exp.Filter{iI(si)}).ex  ...
%     ) )./( SP.(MD.exp.FPhore{iI(si)}).ab .* l_path + eps );
% 
% 
% 
% c = 0.0002;
% I = @(c)...
%     Gain.*Ceff.*MD.exp.Exposure(iI(si)).*10^-3.*exPower.(MD.exp.Filter{iI(si)})...
%     ./(MD.cam.PixNumX*MD.cam.PixNumY) ...
%     .*SP.(MD.exp.FPhore{iI(si)}).QY.* ...
%     sum( QEcam.* ...
%     SP.(MD.exp.Filter{iI(si)}).em .* ...
%     SP.(MD.exp.FPhore{iI(si)}).em./sum(SP.(MD.exp.FPhore{iI(si)}).em) ) .* ...
%     sum( SP.(MD.exp.Filter{iI(si)}).ex./sum(SP.(MD.exp.Filter{iI(si)}).ex) .* ...
%     (1 - 10.^(-SP.(MD.exp.FPhore{iI(si)}).ab .* ...
%                SP.(MD.exp.FPhore{iI(si)}).mec.* l_path .* c) )  )  ;
% 










end