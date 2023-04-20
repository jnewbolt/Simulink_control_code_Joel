%% Parameter sweep

% cd('R:\ENG_Breuer_Shared\ehandyca\handy_simulink_control_code')

if ~exist('experiment','var') ||  ~exist('bias_unloaded','var') || ~exist('bias_loaded','var')
    error('Run "setup_DAQ_simulink" to establish experimental setup. Vars "experiment", "bias_unloaded", "bias_loaded" must be established.')
end

% Date (don't auto-generate date in case experiment runs overnight)
date_start = '20230420';

% Save folder location
% FOLDERNAME = (['R:\ENG_Breuer_Shared\ehandyca\DATA_main_repo\',date_start,'_TandemTuesday_4c_separation_3alphaSweep_diffAlphaValues_APHPH_A3E_02']);
FOLDERNAME = (['D:\Experiments\2foil\',date_start]);

mkdir(FOLDERNAME);

%% Sweep parameters

% non-changing parameters
U = 0.33;
phi = -90;
num_cyc = 20;
transient_cycs = 3;
fred = 0.11;
freq = fred*U/foil.chord;
% freq = 0.65; % very close ~0.649

% non-dim parameters
P1star_vec = [50]; %,60,80];
H1star = 1.2;
P2star_vec = 50; %70; % 65,75
H2star_vec = 0.6; %[0.6,0.8,1.0,1.2,1.4,1.6];
% H2star_vec = [0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2];
% H2star_vec = [1.4,1.6,1.8,2.0,2.2];

phase_vec = 0; %[-180,-120,-60,0,60,120,180];

%% Experimental loop

for P1star = P1star_vec
    for P2star = P2star_vec
        for H2star = H2star_vec
            for phase = phase_vec

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
                heave1 = H1star*foil.chord;
                pitch2 = P2star;
                heave2 = H2star*foil.chord;
                
                freq = fred*U/foil.chord;
                freq = 0.65;

                aT4 = atan(-2*pi*(heave1/foil.chord)*fred) + deg2rad(pitch1);

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
                % Assemble output profiles
                profs(:,1) = [zeros(experiment.motion_delay,1); ramp_p1; pprof1+experiment.offset_home(1)+bias_trial.pitch(1); flip(ramp_p1)];
                profs(:,2) = [zeros(experiment.motion_delay,1); ramp_h1; pprof2+experiment.offset_home(2)+bias_trial.heave(1); flip(ramp_h1)];
                profs(:,3) = [ramp_p2; pprof3+experiment.offset_home(3)+bias_trial.pitch(2); flip(ramp_p2); zeros(experiment.motion_delay,1)];
                profs(:,4) = [ramp_h2; pprof4+experiment.offset_home(4)+bias_trial.heave(2); flip(ramp_h2); zeros(experiment.motion_delay,1)];
                profs(:,5) = [zeros(experiment.motion_delay,1); zeros(size(ramp_p1)); rprof5; zeros(size(ramp_p1))]; % reference signal
                
                % plot trajectories
                plot_profiles(profs);
                
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
                
                % clear variables before next experiment
                clear raw_encoder_p1 raw_encoder_h1 raw_encoder_p2 raw_encoder_h2 raw_force_wallace raw_force_gromit ref_signal

                set_param('simulink_traverse_control','SimulationCommand','start');
                
                disp('Running traverse...')
                pause(sim_time);
                disp('Acquiring data...')
                while ~exist('raw_encoder_p1','var') || ~exist('raw_encoder_h1','var') || ~exist('raw_encoder_p2','var') || ~exist('raw_encoder_h2','var') || ~exist('raw_force_wallace','var') || ~exist('raw_force_gromit','var') || ~exist('ref_signal','var')
                        pause(5)
                        disp('Loading...')
                end
                disp('Done')
                                
                %% Data conversion
                
                range_x = find(ref_signal);
                raw_encoders = [raw_encoder_p1, raw_encoder_h1, raw_encoder_p2, raw_encoder_h2];
                out = convert_output(raw_encoders, raw_force_wallace, raw_force_gromit, raw_vectrino, ref_signal, bias_trial, range_x, experiment.offset_home);
                
                %% Save data

                motor_warning_flag = 0;
%                 FILENAME = (['\',date_start,'_TandemTuesday_4c_separation_3alphaSweep_diffAlpha_',...
%                     'aT4=',num2str(aT4,3),'_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);
                FILENAME = (['\',date_start,'_PrescribedMotion_',...
                    'aT4=',num2str(aT4,3),'_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);

                save(fullfile(FOLDERNAME,FILENAME));

                %% Check for misalignment
                
                if abs(mean(out(:,18))) > 0.08 % if average value of symmetric signal is more than 0.5 deg
                    warning('Gromit pitch motor (Hudson) jerked. Foil will be realigned and trial will be repeated');
                    motor_warning_flag = 1; % raises flag if misalignment due to jerk was detected
%                     FILENAME = ([date_start,'_TandemTuesday_4c_separation_3alphaSweep_diffAlpha_',...
%                         'aT4=',num2str(aT4,3),'_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);
                    FILENAME = (['\',date_start,'_PrescribedMotion_',...
                    'aT4=',num2str(aT4,3),'_p2=',num2str(pitch2,2),'deg_h2=',num2str(heave2/foil.chord,3),'c_ph=',num2str(phase),'deg.mat']);

                    % realign gromit
                    traverse = 'g';
                    run("find_zero_pitch_simulink.m")

                end
                close all
            end
        end
    end
end

% message = ['The experiment finished at ',string(datetime),'. Come and check it out!'];
% sendmail('joel_newbolt@brown.edu','Experiment done',message);

disp('End of experiment')
