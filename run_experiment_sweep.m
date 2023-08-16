%% Parameter sweep
% Run this script after "setup_DAQ_simulink.m" in order to run an experiment with trials that change the
% flapping parameters

% check if necessary variables were setup
if ~exist('Parameters','var') ||  ~exist('Biases','var') || ~exist('BiasesLoaded','var')
    error('Run "setup_DAQ_simulink" to establish experimental setup. Vars "experiment", "bias_unloaded", "bias_loaded" must be established.')
end

%% Sweep parameters
% non-changing parameters
flowSpeedMetersPerSec = 0.30;
phaseLagPitchDeg = 90;
nCycles = 30;
nTransientCycs = 3;
% variable parameters
freq = 0.8889; % Frequency in cycles/sec
freqG = freq; freqW = freq;
pitchAmpDegGvec = 0;%10; %(0:10:30); %,60,80]; % pitch amplitude in degrees
heaveAmpMetersGvec = 0;%0.5; %(0:0.2:0.6); % heave amplitude in chord lengths
pitchAmpDegWvec = 0; %70; % 65,75
heaveAmpMetersWvec = 0;%(0:0.05:1.1); %[0.6,0.8,1.0,1.2,1.4,1.6];
initialPhaseLagWbehindG = -180; 
phaseLagStep = 20; % phase change between trials
phaseLagWbehindGvec = -180;%(initialPhaseLagWbehindG:phaseLagStep:180);

nTrials = length(pitchAmpDegGvec)*length(heaveAmpMetersGvec)*length(pitchAmpDegWvec)*length(heaveAmpMetersWvec)*length(phaseLagWbehindGvec);
iTrial = 1;
%% Move motors to starting position
disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure they have clearance then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,Parameters.firstFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,Parameters.firstFoilHeaveOffsetMeters); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(0,Parameters.secondFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(0,Parameters.secondFoilHeaveOffsetMeters); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Experimental loop
for pitchAmpDegG = pitchAmpDegGvec
for heaveAmpMetersG = heaveAmpMetersGvec
for pitchAmpDegW = pitchAmpDegWvec
for heaveAmpMetersW = heaveAmpMetersWvec
    phaseLagStep = initialPhaseLagWbehindG; % initial value of phase difference in degrees
    while phaseLagStep <= max(phaseLagWbehindGvec) % while loop to change the phase so trial can be repeated if the simulink model fails to run
        %% Take loaded bias measurement
        BiasesLoaded = Biases; % use bias_loaded as the bias to update the pitch and heave biases in the "fin_bias_simulink" routine
        run("find_bias_simulink.m"); % find another loaded bias that contains the drifting and the load bias
        BiasesNewLoaded = Biases; % establish new loaded bias
        % for experiment: bias_trial = bias_newloaded - bias_loaded + bias_unloaded
        BiasesTrial.Wallace = BiasesNewLoaded.Wallace - BiasesLoaded.Wallace + BiasesNoLoad.Wallace;
        BiasesTrial.Gromit = BiasesNewLoaded.Gromit - BiasesLoaded.Gromit + BiasesNoLoad.Gromit;

        %% Check velocity and acceleration limits
        heaveVelocityLimit = 0.5; % Heave velocity limit in m/s
        heaveAccelerationLimit = 4.5; % Heave acceleration limit in m/s/s
        heaveVelocityMaxG = heaveAmpMetersG*2*pi*freqG;
        heaveVelocityMaxW = heaveAmpMetersW*2*pi*freqW;
        heaveAccelerationMaxG = heaveAmpMetersG*(2*pi*freqG)^2;
        heaveAccelerationMaxW = heaveAmpMetersW*(2*pi*freqW)^2;
        if heaveVelocityMaxG > heaveVelocityLimit || heaveVelocityMaxW > heaveVelocityLimit || heaveAccelerationMaxG > heaveAccelerationLimit || heaveAccelerationMaxW > heaveAccelerationLimit
            break % Jump to next trial if command limits are exceeded
        end

        %% Profile generation
        % Generate experiment profiles
        [~, pprof1] = trajectory_experiment(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, pitchAmpDegG, phaseLagPitchDeg, 0);
        [~, pprof2] = trajectory_experiment(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, heaveAmpMetersG, 0, 0);
        [~, pprof3] = trajectory_experiment(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, pitchAmpDegW, phaseLagStep+phaseLagPitchDeg, 0);
        [~, pprof4] = trajectory_experiment(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, heaveAmpMetersW , phaseLagStep, 0);
        [~, rprof5] = trajectory_experiment(nCycles-4, freq, EP.sampleRate, nTransientCycs+2, nTransientCycs+2, 1, 0, 1); % reference signal
        
        clear trajPitchDegreesG trajHeaveMetersG trajPitchDegreesW trajHeaveMetersW trajRefSig
        % Assemble output profiles
        trajPitchDegreesG = endPitchDegG+[zeros(EP.heaveDelayW+EP.pitchDelayW,1); pprof1];
        trajHeaveMetersG = endHeaveMetersG+[zeros(EP.heaveDelayW+EP.pitchDelayW,1); pprof2];
        trajPitchDegreesW = endPitchDegW+[zeros(EP.heaveDelayW,1); pprof3; zeros(EP.pitchDelayW,1)];
        trajHeaveMetersW = endHeaveMetersW+[pprof4; zeros(EP.heaveDelayW+EP.pitchDelayW,1)];
        trajRefSig = [zeros(EP.heaveDelayW,1); zeros(EP.pitchDelayW,1); rprof5]; % reference signal
        
        % plot trajectories
        plot_profiles(trajPitchDegreesG,trajHeaveMetersG,trajPitchDegreesW,trajHeaveMetersW);
        
        % convert into time series to be output to simulink
        % convert into time series to be output to simulink
        times = (0:length(trajPitchDegreesG)-1)'/EP.sampleRate; % time vector to create time series objects
        pitchDegreesG = timeseries(trajPitchDegreesG,times);
        heaveMetersG = timeseries(trajHeaveMetersG,times);
        pitchDegreesW = timeseries(trajPitchDegreesW,times);
        heaveMetersW = timeseries(trajHeaveMetersW,times);
        syncSig = timeseries(trajRefSig,times);

        % experiment simulation time
        simTime = ceil(times(end))+2;
        disp(['Expected simulation time: ', num2str(simTime), ' seconds']);
        
        %% Run actual experiment
        simStatus = 'stopped';
        while strcmp(simStatus,'stopped')
            % clear variables before next experiment
            clear rawEncoderPitchCountsG rawEncoderHeaveCountsG rawEncoderPitchCountsW rawEncoderHeaveCountsW ...
            rawForceVoltsW rawForceVoltsG   rawVoltsVectrino rawVoltsAccelmeter refSig

            disp(['Beginning trial ',num2str(iTrial),' of ',num2str(nTrials)])
            set_param('simulink_traverse_control','SimulationCommand','start');
            simStatus = get_param('simulink_traverse_control','SimulationStatus');
            disp('Running traverse...')
            pause(simTime+5);
%                     disp('Acquiring data...')
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
        
        disp('Done')
                        
        %% Data conversion
        
        rangeTimes = find(refSig);
        rawEncoders = [rawEncoderPitchCountsG, rawEncoderHeaveCountsG, rawEncoderPitchCountsW, rawEncoderHeaveCountsW];
        Measurements = convert_output(rawEncoders, rawForceVoltsW, rawForceVoltsG, rawVoltsVectrino, rawVoltsAccelmeter, refSig, Biases, rangeTimes, EP);

        %% Save data
        trialfilename = ['trial_',num2str(iTrial)];
        save([folder_name,'\',trialfilename]);

        %% Check for misalignment
        
%                 if abs(mean(out(:,18))) > 0.08 % if average value of symmetric signal is more than 0.5 deg
%                     warning('Gromit pitch motor (Hudson) jerked. Foil will be realigned and trial will be repeated');
%                     motor_warning_flag = 1; % raises flag if misalignment due to jerk was detected
% %                     FILENAME = ([date_start,'_TandemTuesday_4c_separation_3alphaSweep_diffAlpha_',...
% %                         'aT4=',num2str(aT4,3),'_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);
%                     FILENAME = (['\',date_start,'_PrescribedMotion',...
%                     '_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);
% 
%                     % realign gromit
%                     traverse = 'g';
%                     run("find_zero_pitch_simulink.m")
% 
%                 end
%                 close all
%% Step forward the phase
    phaseLagStep = phaseLagStep + phaseLagStep;
    iTrial = iTrial+1;
    end
end
end
end
end

%% Move motors to home position
% check for required parameters:
if ~exist('Parameters','var')
    error('Missing necessary variables from workspace')
end

disp('The traverses will move to their home positions.')
% Run the move to start positions
[startPitchDegG, endPitchDegG] = deal(endPitchDegG,0); 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0);
[startPitchDegW, endPitchDegW] = deal(endPitchDegW,0); 
[startHeaveMetersW, endHeaveMetersW] = deal(endHeaveMetersW,0);
run('move_to_position')
clearvars -except Parameters Measurements Biases

%% Send email
% send an email letting me know the experiment was completed
message = strjoin(['The experiment finished at ',string(datetime),'. Come and check it out!']);
sendmail('joel_newbolt@brown.edu','Experiment done',message);

disp('End of experiment')
