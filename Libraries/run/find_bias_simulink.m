% function [Measurements, Biases] = find_bias_simulink(EP)
%% Find bias - simulink control
    % Find bias in the force measurements
    % NOTE: only overwrites the force biases unless bias is not specified
    % bias is an optional argument
EP = ExperimentParameters; 
sampleTime = 1/EP.sampleRate;
%% Stationary command profiles
% ramp to offset time
rampTime = 5; % in [s]
% Generate ramp profiles
[~, rampPitchG, rampHeaveG] = ramp_fn(rampTime, EP, 'Gromit');
[~, rampPitchW, rampHeaveW] = ramp_fn(rampTime, EP, 'Wallace');

% Zero motion profile
trajDuration = 20; % length of time of bias measurement in seconds
trajStationary = zeros(trajDuration*EP.sampleRate,1);
motionDelayTimesteps = EP.motionDelay*(1000/EP.sampleRate);

clear profs % in case this scripts is run after a previous one
% gromit pitch
trajPitchDegreesG = [zeros(motionDelayTimesteps,1); rampPitchG; trajStationary+EP.firstFoilPitchOffsetDegrees; flip(rampPitchG)];
% gromit heave
trajHeaveMetersG = [zeros(motionDelayTimesteps,1); rampHeaveG; trajStationary+EP.firstFoilHeaveOffsetMeters; flip(rampHeaveG)];
% wallace pitch
trajPitchDegreesW = [rampPitchW; trajStationary+EP.secondFoilPitchOffsetDegrees; flip(rampPitchW); zeros(motionDelayTimesteps,1)];
% wallace heave
trajHeaveMetersW = [rampHeaveW; trajStationary+EP.secondFoilHeaveOffsetMeters; flip(rampHeaveW); zeros(motionDelayTimesteps,1)];

% reference signal - has 0 for ramps and delays and 1 for usable data
trajRefSig = [zeros(motionDelayTimesteps,1); zeros(size(rampPitchG)); ones(size(trajStationary)); zeros(size(rampPitchG))];
% NOTE: the profile for sync_sig should be value 1 some time AFTER the ramp is finished,
% and should go back to 0 some time BEFORE the 0 profile ends.

% plot trajectories
%     plot_profiles(profs);

% convert into time series to be output to simulink
times = (0:length(trajPitchDegreesG)-1)'/EP.sampleRate; % time vector to create time series objects
pitchDegreesG = timeseries(trajPitchDegreesG,times);
heaveMetersG = timeseries(trajHeaveMetersG,times);
pitchDegreesW = timeseries(trajPitchDegreesW,times);
heaveMetersW = timeseries(trajHeaveMetersW,times);
syncSig = timeseries(trajRefSig,times);

% simulation time
simTime = ceil(times(end))+2;
disp(['Expected simulation time: ', num2str(simTime), ' seconds']);

% pass parameters for gromit heave gain in simulation
freq = 0; heave1 = 0;

    %% Run traverse
simStatus='stopped'; % Check if the model ran correctly
while strcmp(simStatus,'stopped') 
    % clear variables in case done after an experiment
%         clear raw_encoder_p1 raw_encoder_h1 raw_encoder_p2 raw_encoder_h2 raw_force_wallace raw_force_gromit ref_signal
        clear rawEncoderPitchCountsG rawEncoderHeaveCountsG rawEncoderPitchCountsW rawEncoderHeaveCountsW ...
            rawForceVoltsW rawForceVoltsG   rawVoltsVectrino rawVoltsAccelmeter refSig
   
%     set_param('simulink_traverse_control','sampleTime',);
    set_param('simulink_traverse_control','SimulationCommand','start');
    simStatus = get_param('simulink_traverse_control','SimulationStatus');
    disp('Running traverse...')
    pause(simTime+5);
%     disp('Acquiring data...')
    while ~exist('refSig','var')
            simStatus = get_param('simulink_traverse_control','SimulationStatus');
            pause(5)
            disp('Loading...')
            if strcmp(simStatus,'stopped') % Break out if the trial is taking unexpectedly long to avoid hanging on failure to build model
                delete('SimulationCache\simulink_traverse_control.slxc')
                disp('Model failed to run.  Deleting cache (.slxc) and trying again.')
                break
            end
    end
end
%     sim_status = 'running';
%     while ~strcmp(sim_status,'stopped')
%         sim_status = get_param('simulink_traverse_control','SimulationStatus');
%         pause(5);
%         disp('Loading...')
%     end
% 
%     if exist('raw_encoder_p1','var') && exist('raw_encoder_h1','var') && exist('raw_encoder_p2','var') && exist('raw_encoder_h2','var') && exist('raw_force_wallace','var') && exist('raw_force_gromit','var') && exist('ref_signal','var')
%         disp('All data acquired')
%     end

disp('Done')

%% Calculate force biases

rangeTimes = find(refSig);
rangeTimes = rangeTimes(1000:end);

Biases.forceVoltsW = mean(rawForceVoltsW(rangeTimes,:),1);
Biases.forceVoltsStdevW = std(rawForceVoltsW(rangeTimes,:),1);

Biases.forceVoltsG = mean(rawForceVoltsG(rangeTimes,:),1);
Biases.forceVoltsStdevG = std(rawForceVoltsG(rangeTimes,:),1);

Biases.accmeterVolts = mean(rawVoltsAccelmeter(rangeTimes,:),1);
Biases.accmeterVoltsStdev = std(rawVoltsAccelmeter(rangeTimes,:),1);

    %% Convert output

rawEncoders = [rawEncoderPitchCountsG, rawEncoderHeaveCountsG, rawEncoderPitchCountsW, rawEncoderHeaveCountsW];
Measurements = convert_output(rawEncoders, rawForceVoltsW, rawForceVoltsG, rawVoltsVectrino, rawVoltsAccelmeter, refSig, Biases, rangeTimes, EP);

%     %% Plot results
% 
%     figure()
%     force_resolution = [1/32 1/32 1/16];
%     torque_resolution = 1/528;
% 
%     subplot(2,1,2)
%     plot(out(20:end-20,7:9)./force_resolution,'.')
%     hold on 
%     plot(out(20:end-20,10:12)/torque_resolution,'.')
%     hold off
%     title('Wallace (last)')
%     ylabel('Forces and Torques (normalized by resolution)')
%     legend('Fy','Fx','Fz','Ty','Tx','Tz')
% 
%     subplot(2,1,1)
%     plot(out(20:end-20,17:19)./force_resolution,'.')
%     hold on 
%     plot(out(20:end-20,20:22)/torque_resolution,'.')
%     hold off
%     title('Gromit (middle)')
%     ylabel('Forces and Torques (normalized by resolution)')
%     legend('Fy','Fx','Fz','Ty','Tx','Tz')
%     
%% Raise warnings for potentially faulty measurements

if Biases.accmeterVolts < 1.5 || Biases.accmeterVolts > 1.8
    disp(['Accelerometer voltage is ', num2str(Biases.accmeterVolts,3), ' Volts'])
    disp('Warning: Accelerometer voltage outside expected range. Try power cycling.')
end

% Biases.forceVoltsRMSEW = sqrt(mean((rawForceVoltsW- repmat(Biases.forceVoltsW,numel(rawForceVoltsW(:,1)),1)).^2));
% 
% Biases.forceVoltsRMSEG = sqrt(mean((rawForceVoltsG - repmat(Biases.forceVoltsG,numel(rawForceVoltsG(:,1)),1)).^2));
% if sum(bias.RMSEW>[.15 .15 .3 .1 .1 .1])>0 
%         disp(bias.RMSEW);
%     disp('Warning: Wallace error signal above normal. Check wiring/ grounding.')
% end
% if sum(bias.RMSEG>[.15 .15 .3 .1 .1 .1])>0 
%     disp('Warning: Gromit error signal above normal. Check wiring/ grounding.')
% end

% percentFullrangeErrorW = Biases.forceVoltsRMSEW./[660 660 1980 60 60 60]*100;
% percentFullrangeErrorG = Biases.forceVoltsRMSEG./[660 660 1980 60 60 60]*100;

%% Save bias into a subfolder within the active saving folder

time = clock;
dataFolderName = [EP.dataFolderName,'\data'];

% Check how many bias files have been created in order to name the current one being generated
biasFiles = dir([dataFolderName,'\Biases*']);
numBiasFiles = numel(biasFiles);
if numBiasFiles == 0
    filename = [dataFolderName,'\BiasesNoLoad'];
    BiasesNoLoad = Biases;
else
    filename = [dataFolderName,'\BiasesLoaded_',num2str(numBiasFiles)];
    BiasesLoaded = Biases;
end
save(filename);

% end