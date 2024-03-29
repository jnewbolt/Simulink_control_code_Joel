%% Parameter sweep
% Run this script after "setup_DAQ_simulink.m" in order to run an experiment with trials that change the
% flapping parameters

% check if necessary variables were setup
if ~exist('experiment','var') ||  ~exist('bias_unloaded','var') || ~exist('bias_loaded','var')
    error('Run "setup_DAQ_simulink" to establish experimental setup. Vars "experiment", "bias_unloaded", "bias_loaded" must be established.')
end

%% Sweep parameters

% non-changing parameters
U = 0.30;
phi = -90;
num_cyc = 30;
transient_cycs = 3;
% fred = 0.11;
% freq = fred*U/foil.chord;
% freq = 0.65; % very close ~0.649

% non-dim parameters
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
%% Experimental loop

for P1star = P1star_vec
for H1star = H1star_vec
    for P2star = P2star_vec
        for H2star = H2star_vec
            phase = initial_phase; % initial value of phase difference in degrees
            while phase <= max(phase_vec) % while loop to change the phase so trial can be repeated if the simulink model fails to run
                %% Take experiment bias measurement
                
                bias = bias_loaded; % use bias_loaded as the bias to update the pitch and heave biases in the "fin_bias_simulink" routine
                run("find_bias_simulink.m"); % find another loaded bias that contains the drifting and the load bias
                bias_newloaded = bias; % establish new loaded bias
                % for experiment: bias_trial = bias_newloaded - bias_loaded + bias_unloaded
                bias_trial = bias; % inherits pitch and heave biases
                bias_trial.Wallace = bias_newloaded.Wallace - bias_loaded.Wallace + bias_unloaded.Wallace;
                bias_trial.Gromit = bias_newloaded.Gromit - bias_loaded.Gromit + bias_unloaded.Gromit;


                %% Dimensional parameters

                pitch1 = P1star;
%                 heave1 = H1star*foil.chord;
                heave1 = H1star*chord_foil;%0.024*2.0362; % heave of upstream foil in meters
                pitch2 = P2star;
                heave2 = H2star*0.0238; % manual value of cross-stream thickness D for ELLIPTICAL CYLINDER
                
%                 freq = fred*U/foil.chord;
                freq = 0.8889; % Frequency in cycles/sec
                freq1 = freq;
                freq2 = freq;

%                 aT4 = atan(-2*pi*(heave1/foil.chord)*fred) + deg2rad(pitch1);

                %% Check velocity and acceleration limits
                heavevel_limit = 0.5; % Heave velocity limit in m/s
                heaveacc_limit = 4.5; % Heave acceleration limit in m/s/s
                heavevel1_max = heave1*2*pi*freq1;
                heavevel2_max = heave2*2*pi*freq2;
                heaveacc1_max = heave1*(2*pi*freq1)^2;
                heaveacc2_max = heave2*(2*pi*freq2)^2;
                if heavevel1_max > heavevel_limit || heavevel2_max > heavevel_limit || heaveacc1_max > heaveacc_limit || heaveacc2_max > heaveacc_limit
                    break
                end

                %% Profile generation
                
                % ramp to offset time
                ramp_time = 5; % in [s]
                
                % Generate ramp profiles
                [~, ramp_p1, ramp_h1] = ramp_fn(ramp_time, experiment.T, bias_trial, experiment.offset_home, 'g');
                [~, ramp_p2, ramp_h2] = ramp_fn(ramp_time, experiment.T, bias_trial, experiment.offset_home, 'w');
                
                % Generate experiment profiles
                [~, pprof1] = generate_profile(num_cyc, freq, experiment.srate, transient_cycs, transient_cycs, pitch1, phi, 0);
                [~, pprof2] = generate_profile(num_cyc, freq, experiment.srate, transient_cycs, transient_cycs, heave1, 0, 0);
                [~, pprof3] = generate_profile(num_cyc, freq, experiment.srate, transient_cycs, transient_cycs, pitch2, phase+phi, 0);
                [~, pprof4] = generate_profile(num_cyc, freq, experiment.srate, transient_cycs, transient_cycs, heave2, phase, 0);
                [~, rprof5] = generate_profile(num_cyc-4, freq, experiment.srate, transient_cycs+2, transient_cycs+2, 1, 0, 1); % reference signal
                
                clear profs
                heaveW_delay = 50; % phase difference between pitch wallace and heave wallace (heave lags behind pitch)
                % Assemble output profiles
                profs(:,1) = [zeros(heaveW_delay,1); zeros(experiment.motion_delay,1); ramp_p1; pprof1+experiment.offset_home(1)+bias_trial.pitch(1); flip(ramp_p1)];
                profs(:,2) = [zeros(heaveW_delay,1); zeros(experiment.motion_delay,1); ramp_h1; pprof2+experiment.offset_home(2)+bias_trial.heave(1); flip(ramp_h1)];
                profs(:,3) = [zeros(heaveW_delay,1); ramp_p2; pprof3+experiment.offset_home(3)+bias_trial.pitch(2); flip(ramp_p2); zeros(experiment.motion_delay,1)];
                profs(:,4) = [ramp_h2; pprof4+experiment.offset_home(4)+bias_trial.heave(2); flip(ramp_h2); zeros(experiment.motion_delay,1); zeros(heaveW_delay,1)];
                profs(:,5) = [zeros(heaveW_delay,1); zeros(experiment.motion_delay,1); zeros(size(ramp_p1)); rprof5; zeros(size(ramp_p1))]; % reference signal
                
%                 % plot trajectories
%                 plot_profiles(profs);
                
                % convert into time series to be output to simulink
                toime = (0:size(profs,1)-1)'/experiment.srate; % time vector to create time series objects
                outp1 = timeseries(profs(:,1),toime);
                outh1 = timeseries(profs(:,2),toime);
                outp2 = timeseries(profs(:,3),toime);
                outh2 = timeseries(profs(:,4),toime);
                sync_sig = timeseries(profs(:,5),toime);

                % experiment simulation time
                sim_time = ceil(toime(end))+2;
                disp(['Expected simulation time: ', num2str(sim_time), ' seconds']);
                
                %% Run actual experiment
                simStatus = 'stopped';
                while strcmp(simStatus,'stopped')
                    % clear variables before next experiment
                    clear raw_encoder_p1 raw_encoder_h1 raw_encoder_p2 raw_encoder_h2 raw_force_wallace raw_force_gromit ref_signal
                    
                    disp(['Beginning trial ',num2str(trial_number),' of ',num2str(num_trials)])
                    set_param('simulink_traverse_control','SimulationCommand','start');
                    simStatus = get_param('simulink_traverse_control','SimulationStatus');
                    disp('Running traverse...')
                    pause(sim_time+5);
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
                
                range_x = find(ref_signal);
                raw_encoders = [raw_encoder_p1, raw_encoder_h1, raw_encoder_p2, raw_encoder_h2];
                out = convert_output(raw_encoders, raw_force_wallace, raw_force_gromit, raw_vectrino, raw_accelmeter, ref_signal, bias_trial, range_x, experiment.offset_home);
                
                %% Save data

                motor_warning_flag = 0;
%                 FILENAME = (['\',date_start,'_TandemTuesday_4c_separation_3alphaSweep_diffAlpha_',...
%                     'aT4=',num2str(aT4,3),'_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);
%                 FILENAME = (['\',date_start,'_PrescribedMotion',...
%                     '_p1=',num2str(pitch1,2),'deg_h1=',num2str(heave1,3),'cm_ph=',num2str(phase),'deg.mat']);
%                 FILENAME = (['\',date_start,'_PrescribedMotion',...
%                     '_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2,3),'cm_ph=',num2str(phase),'deg.mat']);
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
%%          Step forward the phase
            phase = phase + phase_step;
            trial_number = trial_number+1;
            end
        end
    end
end
end

% send an email letting me know the experiment was completed
message = strjoin(['The experiment finished at ',string(datetime),'. Come and check it out!']);
sendmail('joel_newbolt@brown.edu','Experiment done',message);

disp('End of experiment')
