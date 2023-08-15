function M = convert_output(rawEncoders, rawForceVoltsW, rawForceVoltsG, rawVoltsVectrino, rawVoltsAccelmeter, ...
    refSig, Biases, rangeTimes, EP)
% % Required inputs:
% raw_encoders
% raw_wallace
% raw_gromit
% % Optional inputs:
% bias
% range_x
% offset_p2
% offset_h2

%% Encoders

if ~exist('rangeTimes', 'var') || isempty(rangeTimes)
    rangeTimes = 1:size(rawEncoders,1);
end

if ~exist('Biases','var') || isempty(Biases)
    Biases.forceVoltsW = 0;
    Biases.forceVoltsG = 0;
end

if ~exist('offset_home','var') || isempty(offset_home)
    offset_p2 = 0;
    offset_h2 = 0;
else
    offset_p2 = offset_home(3);
    offset_h2 = offset_home(4);
end

rawEncoders(rawEncoders(:,1:2)>1e6) = rawEncoders(rawEncoders(:,1:2)>1e6)-2^32;

M.pitchDegreesG = rawEncoders(rangeTimes,1).*360/10000;
M.heaveMetersG = -rawEncoders(rangeTimes,2).*0.0254/8000;
M.pitchDegreesW = conv_encodertheta_traverse(rawEncoders(rangeTimes,3), 0);
M.heaveMetersW = conv_encodery_traverse(rawEncoders(rangeTimes,4), 0);

% refSig = ref_signal(rangeTimes);

%% Force transducers

M.forcesNewtonsW = conv_output_wallace_force(rawForceVoltsW(rangeTimes,:),Biases.forceVoltsW);
M.forcesNewtonsG = conv_output_gromit_force(rawForceVoltsG(rangeTimes,:),Biases.forceVoltsG);

%% Vectrino

M.flowMetersPerSecondVectrino = (rawVoltsVectrino(rangeTimes,:)-2.5)*2/5;

%% Accelerometer

forceSensorMassKilograms = 0.6; % Mass of mount 
accmeterMetersPerSecondSqPerVolt = 9.81; % Meters per seconds squared per volt
M.forceInertialLoadG = (forceSensorMassKilograms+EP.Foils.firstFoil.mass)*accmeterMetersPerSecondSqPerVolt* ...
    (rawVoltsAccelmeter(rangeTimes,:)-Biases.accmeterVolts); 

end

