%% Parameter sweep
% Run this script after "setup_DAQ_simulink.m" in order to run an experiment with trials that change the
% flapping parameters
tic % begin timing the experiment
% check if necessary variables were setup
if ~exist('Parameters','var') ||  ~exist('Biases','var')
    error('Run "setup_DAQ_simulink" to establish experimental setup. Structures "Parameters" and "Biases" must be established.')
end
%% Move motors to starting position
disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure MOTORS have CLEARANCE and FLUME IS ON, then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,Parameters.pitchOffsetDegG); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,Parameters.heaveOffsetMetersG); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(0,Parameters.pitchOffsetDegW); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(0,Parameters.heaveOffsetMetersW); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Biases endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Sweep parameters
% non-changing parameters
flowSpeedMetersPerSec = 0.30;
phaseLagPitchDeg = 90;
nCyclesG = 30;
nTransientCycsG = 2;
% variable parameters
freqW = 0.889;
freqGvec = freqW; 
pitchAmpDegGvec = 0;%10; %(0:10:30); %,60,80]; % pitch amplitude in degrees
heaveAmpMetersGvec = 0;%0.05; %(0:0.2:0.6); % heave amplitude in chord lengths
pitchAmpDegWvec = (0:1:3);%10; %70; % 65,75
heaveAmpMetersWvec = (0:0.01:0.03);%0.05;%(0:0.05:1.1); %[0.6,0.8,1.0,1.2,1.4,1.6];
initialPhaseLagWbehindG = 0; 
phaseLagStep = 20; % phase change between trials
phaseLagWbehindGvec = 0;%(initialPhaseLagWbehindG:phaseLagStep:180);

nTrials = length(pitchAmpDegGvec)*length(heaveAmpMetersGvec)*length(pitchAmpDegWvec)*length(heaveAmpMetersWvec)*length(phaseLagWbehindGvec);
iTrial = 1;

%% Experimental loop
for freqG = freqGvec
for pitchAmpDegG = pitchAmpDegGvec
for heaveAmpMetersG = heaveAmpMetersGvec
for pitchAmpDegW = pitchAmpDegWvec
for heaveAmpMetersW = heaveAmpMetersWvec
    phaseLagWbehindG = initialPhaseLagWbehindG; % initial value of phase difference in degrees
    while phaseLagWbehindG <= max(phaseLagWbehindGvec) % while loop to change the phase so trial can be repeated if the simulink model fails to run
        %% Take loaded bias measurement
        disp('Finding LOADED BIAS')
        run("find_bias_simulink.m"); % find another loaded bias that contains the drifting and the loaded bias
%         clearvars -except Parameters Biases endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW
        % for experiment: BiasesNoLoad = BiasesLoaded - BiasesLoaded0 + BiasesNoLoad0
        Biases.NoLoad = Biases.NoLoad0; % Passes accelerometer bias
        Biases.NoLoad.forceVoltsW = Biases.Loaded.forceVoltsW - Biases.Loaded0.forceVoltsW + Biases.NoLoad0.forceVoltsW;
        Biases.NoLoad.forceVoltsG = Biases.Loaded.forceVoltsG - Biases.Loaded0.forceVoltsG + Biases.NoLoad0.forceVoltsG;    
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
        % Dependent parameters
        nCyclesW = (freqW/freqG)*nCyclesG;
        nTransientCycsW = (freqW/freqG)*nTransientCycsG;
        % Generate experiment profiles
        [~, pprof1] = trajectory_experiment(nCyclesG, freqG, Parameters.sampleRate, nTransientCycsG, nTransientCycsG, pitchAmpDegG, phaseLagPitchDeg, 0);
        [~, pprof2] = trajectory_experiment(nCyclesG, freqG, Parameters.sampleRate, nTransientCycsG, nTransientCycsG, heaveAmpMetersG, 0, 0);
        [~, pprof3] = trajectory_experiment(nCyclesW, freqW, Parameters.sampleRate, nTransientCycsW, nTransientCycsW, pitchAmpDegW, phaseLagWbehindG+phaseLagPitchDeg, 0);
        [~, pprof4] = trajectory_experiment(nCyclesW, freqW, Parameters.sampleRate, nTransientCycsW, nTransientCycsW, heaveAmpMetersW , phaseLagWbehindG, 0);
        [~, rprof5] = trajectory_experiment(nCyclesG, freqG, Parameters.sampleRate, nTransientCycsG, nTransientCycsG, 1, 0, 1); % reference signal
        
        clear trajPitchDegreesG trajHeaveMetersG trajPitchDegreesW trajHeaveMetersW trajRefSig
        % Assemble output profiles
        trajPitchDegreesG = endPitchDegG+[zeros(Parameters.heaveDelayW+Parameters.pitchDelayW,1); pprof1];
        trajHeaveMetersG = endHeaveMetersG+[zeros(Parameters.heaveDelayW+Parameters.pitchDelayW,1); pprof2];
        trajPitchDegreesW = endPitchDegW+[zeros(Parameters.heaveDelayW,1); pprof3; zeros(Parameters.pitchDelayW,1)];
        trajHeaveMetersW = endHeaveMetersW+[pprof4; zeros(Parameters.heaveDelayW+Parameters.pitchDelayW,1)];
        trajRefSig = [zeros(Parameters.heaveDelayW,1); zeros(Parameters.pitchDelayW,1); rprof5]; % reference signal
        
        % plot trajectories
        plot_profiles(trajPitchDegreesG,trajHeaveMetersG,trajPitchDegreesW,trajHeaveMetersW);
        
        % convert into time series to be output to simulink
        % convert into time series to be output to simulink
        times = (0:length(trajPitchDegreesG)-1)'/Parameters.sampleRate; % time vector to create time series objects
        pitchDegreesG = timeseries(trajPitchDegreesG,times);
        heaveMetersG = timeseries(trajHeaveMetersG,times);
        pitchDegreesW = timeseries(trajPitchDegreesW,times);
        heaveMetersW = timeseries(trajHeaveMetersW,times);
        syncSig = timeseries(trajRefSig,times);

        % experiment simulation time
        simTime = ceil(times(end))+2;
        disp(['Expected simulation time: ', num2str(simTime), ' seconds']);
        
        %% Run actual experiment
        freqGain = freqG; heaveGain = heaveAmpMetersG; % These values modify the Gromit motor command to get desired plant positions 
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
        Measurements = convert_output(rawEncoders, rawForceVoltsW, rawForceVoltsG, rawVoltsVectrino, rawVoltsAccelmeter, refSig, Biases.NoLoad, rangeTimes, Parameters);

        %% Save data
        trialfilename = ['trial_',num2str(iTrial)];
        dataFolderName = [Parameters.dataFolderName,'\data'];
        save([dataFolderName,'\',trialfilename]);

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
    phaseLagWbehindG = phaseLagWbehindG + phaseLagStep;
    iTrial = iTrial+1;
    end
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
experimentDuration = toc;
message = strjoin(['The experiment finished at ',string(datetime),' with a total runtime of ',num2str(experimentDuration/60,1),' minutes. Come and check it out!']);
sendmail('joel_newbolt@brown.edu','Experiment done',message);

disp('End of experiment')
