%% Parameter sweep
% Run this script after "setup_DAQ_simulink.m" in order to run an experiment with trials that change the
% flapping parameters

% check if necessary variables were setup
if ~exist('ExperimentParameters','var') ||  ~exist('Biases','var') || ~exist('BiasesLoaded','var')
    error('Run "setup_DAQ_simulink" to establish experimental setup. Vars "experiment", "bias_unloaded", "bias_loaded" must be established.')
end

%% Sweep parameters
% non-changing parameters
flowSpeedMetersPerSec = 0.30;
phi = -90;
nCycles = 30;
nTransientCycs = 3;
% fred = 0.11;
% freq = fred*U/foil.chord;
% freq = 0.65; % very close ~0.649

% non-dim parameters
freq = 0.8889; % Frequency in cycles/sec
P1star_vec = 0;%10; %(0:10:30); %,60,80]; % pitch amplitude in degrees
H1star_vec = 0;%0.5; %(0:0.2:0.6); % heave amplitude in chord lengths
chord_foil = 0.06; % chord length of upstream foil in meters
P2star_vec = 0; %70; % 65,75
H2star_vec = 0;%(0:0.05:1.1); %[0.6,0.8,1.0,1.2,1.4,1.6];
initial_phase = -180; 
phase_step = 20; % phase change between trials
phase_vec = -180;%(initial_phase:phase_step:180);

num_trials = length(P1star_vec)*length(H1star_vec)*length(P2star_vec)*length(H2star_vec)*length(phase_vec);
trial_number = 1;
%% Move motors to starting position
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure they have clearance then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,ExperimentParameters.firstFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(0,ExperimentParameters.secondFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(0,ExperimentParameters.secondFoilHeaveOffsetMeters); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except ExperimentParameters

%% Experimental loop
for P1star = P1star_vec
for H1star = H1star_vec
for P2star = P2star_vec
for H2star = H2star_vec
    phase = initial_phase; % initial value of phase difference in degrees
    while phase <= max(phase_vec) % while loop to change the phase so trial can be repeated if the simulink model fails to run
        %% Take loaded bias measurement
        BiasesLoaded = Biases; % use bias_loaded as the bias to update the pitch and heave biases in the "fin_bias_simulink" routine
        run("find_bias_simulink.m"); % find another loaded bias that contains the drifting and the load bias
        BiasesNewLoaded = Biases; % establish new loaded bias
        % for experiment: bias_trial = bias_newloaded - bias_loaded + bias_unloaded
        BiasesTrial.Wallace = BiasesNewLoaded.Wallace - BiasesLoaded.Wallace + BiasesNoLoad.Wallace;
        BiasesTrial.Gromit = BiasesNewLoaded.Gromit - BiasesLoaded.Gromit + BiasesNoLoad.Gromit;
        %% Dimensional parameters
        pitch1 = P1star;
%                 heave1 = H1star*foil.chord;
        heave1 = H1star*chord_foil;%0.024*2.0362; % heave of upstream foil in meters
        pitch2 = P2star;
        heave2 = H2star*0.0238; % manual value of cross-stream thickness D for ELLIPTICAL CYLINDER        
        freq1 = freq;
        freq2 = freq;

        %% Check velocity and acceleration limits
        heavevel_limit = 0.5; % Heave velocity limit in m/s
        heaveacc_limit = 4.5; % Heave acceleration limit in m/s/s
        heavevel1_max = heave1*2*pi*freq1;
        heavevel2_max = heave2*2*pi*freq2;
        heaveacc1_max = heave1*(2*pi*freq1)^2;
        heaveacc2_max = heave2*(2*pi*freq2)^2;
        if heavevel1_max > heavevel_limit || heavevel2_max > heavevel_limit || heaveacc1_max > heaveacc_limit || heaveacc2_max > heaveacc_limit
            break % Jump to next trial if command limits are exceeded
        end

        %% Profile generation
        % Generate experiment profiles
        [~, pprof1] = generate_profile(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, pitch1, phi, 0);
        [~, pprof2] = generate_profile(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, heave1, 0, 0);
        [~, pprof3] = generate_profile(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, pitch2, phase+phi, 0);
        [~, pprof4] = generate_profile(nCycles, freq, EP.sampleRate, nTransientCycs, nTransientCycs, heave2, phase, 0);
        [~, rprof5] = generate_profile(nCycles-4, freq, EP.sampleRate, nTransientCycs+2, nTransientCycs+2, 1, 0, 1); % reference signal
        
        clear profs
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
            clear raw_encoder_p1 raw_encoder_h1 raw_encoder_p2 raw_encoder_h2 raw_force_wallace raw_force_gromit ref_signal
            
            disp(['Beginning trial ',num2str(trial_number),' of ',num2str(num_trials)])
            set_param('simulink_traverse_control','SimulationCommand','start');
            simStatus = get_param('simulink_traverse_control','SimulationStatus');
            disp('Running traverse...')
            pause(simTime+5);
%                     disp('Acquiring data...')
            while ~exist('raw_encoder_p1','var') || ~exist('raw_encoder_h1','var') || ~exist('raw_encoder_p2','var') || ~exist('raw_encoder_h2','var') || ~exist('raw_force_wallace','var') || ~exist('raw_force_gromit','var') || ~exist('ref_signal','var')
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
        trialfilename = ['trial_',num2str(trial_number)];
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
    phase = phase + phase_step;
    trial_number = trial_number+1;
    end
end
end
end
end

%% Move motors to home position
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp('The traverses will move to their home positions.')
% Run the move to start positions
[startPitchDegG, endPitchDegG] = deal(endPitchDegG,0); 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0);
[startPitchDegW, endPitchDegW] = deal(endPitchDegW,0); 
[startHeaveMetersW, endHeaveMetersW] = deal(endHeaveMetersW,0);
run('move_to_position')
clearvars -except ExperimentParameters Measurements Biases

%% Send email
% send an email letting me know the experiment was completed
message = strjoin(['The experiment finished at ',string(datetime),'. Come and check it out!']);
sendmail('joel_newbolt@brown.edu','Experiment done',message);

disp('End of experiment')
