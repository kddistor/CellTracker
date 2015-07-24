%IMAN_FRETFACTOR_EST
%   Estimate Fret Factor

function ff = iman_fretfactor_est(MD, SP)

%Default Fret pair = CFP/YFP
cD = 'cfp';  cA = 'yfp';

%Ensure spectral parameters are available
if ~exist('SP', 'var') || isempty(SP)
    %Run routine to extract spectral parameters
    SP = iq_getspectralpar;
end

%Get excitation power estimates
exPower = iq_excitation_powermap(MD.exp.Filter, MD.exp.ExVolt);

%Get camera QE (quantum efficiency) curve
QEcam = iq_getcameraqe(MD.cam.Desc, SP.WaveLength);


%Get Channel IDs for FRET pair
iD = find( ~cellfun('isempty', regexpi(MD.exp.Channel, cD, 'start')) );
iA = find( ~cellfun('isempty', regexpi(MD.exp.Channel, cA, 'start')) );


%Fret Factor Calculation
% ff = Abs(D_d) * QY(D) * Em(D_d) * t(d) * P(d) ./ ...
%      Abs(A_a) * QY(A) * Em(A_a) * t(a) * P(a);

ff = sum( SP.(MD.exp.Filter{iD}).ex .* SP.(MD.exp.FPhore{iD}).ab ) .* ...
     sum( SP.(MD.exp.Filter{iD}).em .* SP.(MD.exp.FPhore{iD}).em .* QEcam ) .* ...
     SP.(MD.exp.FPhore{iD}).QY .* MD.exp.Exposure(iD) .* exPower.(MD.exp.Filter{iD}) ...
     ./ ...
   ( sum( SP.(MD.exp.Filter{iA}).ex .* SP.(MD.exp.FPhore{iA}).ab ) .* ...
     sum( SP.(MD.exp.Filter{iA}).em .* SP.(MD.exp.FPhore{iA}).em .* QEcam ) .* ...
     SP.(MD.exp.FPhore{iA}).QY .* MD.exp.Exposure(iA) .* ...
     exPower.(MD.exp.Filter{iA}) + eps );  %eps to prevent NaN


end