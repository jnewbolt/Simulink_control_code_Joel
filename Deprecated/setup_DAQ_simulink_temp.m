%% Simulink control code daq parameters

addpath(genpath('Libraries'))


%% DAQ setup parameters

srate = 1000; % sampling frequency
T = 1/srate; % sampling period

motion_delay = 35; % time by which gromit's motion leads the new traverse in [ms]

%% Experimental parameters

f_pump = 17.9; % pump frequency in Hz

U = 0.33; % expected freestream flowspeed
foiltype = 'A3E';
foil = foils_database('A3E');
chord = foil.chord;

num_cyc = 10; % number of full-amplitude experimental cycles
phi = -90; % phase by which pitch motion leads heave motion

fred = 0.12; % reduced frequency
freq = fred*U/chord; % main experimental frequency
freq = 0.65; % won't change the results (0.6492 ~ 0.65)
phase = -60; % inter-foil-phase

pitch1 = 50; % leading pitch amplitude
pitch2 = 80; % trailing pitch amplitude
hstar1 = 0.8; % non-dim leading heave amplitude
hstar2 = 0.8; % non-dim trailing heave amplitude
heave1 = hstar1*chord;
heave2 = hstar2*chord;

%% Control parameters

bias_pitch1 = 0; % this will have to be updated in some other way
bias_heave1 = 0; % this will have to be updated in some other way
bias_pitch2 = 0; % this will have to be updated in some other way
bias_heave2 = 0; % this will have to be updated in some other way

% Ramp function

ramp_time = 5; % time for ramping up to position in seconds\

offset_h2 = 0.25; % offset position for the heave axis of the new traverse
offset_p2 = 180; % offset position for the pitch axis of the new traverse

Npts_ramp = ramp_time/T;
ramp_h2 = offset_h2*(0.5*(1-cos( pi*(0:Npts_ramp-1)/Npts_ramp)))';
ramp_p2 = offset_p2*(0.5*(1-cos( pi*(0:Npts_ramp-1)/Npts_ramp)))';

[t1, pprof1] = generate_profile(num_cyc, freq, srate, 3, 3, pitch1, phi, 0);
[t2, hprof1] = generate_profile(num_cyc, freq, srate, 3, 3, heave1, 0, 0);
[t3, pprof2] = generate_profile(num_cyc, freq, srate, 3, 3, pitch2, phase+phi, 0);
[t4, hprof2] = generate_profile(num_cyc, freq, srate, 3, 3, heave2, phase, 0);

t_ramp = T*(0:Npts_ramp-1); % dimensional ramping time vector

% leading motion profiles
data_p1 = [zeros(length(ramp_p2)+motion_delay,1); pprof1; zeros(length(ramp_p2)-motion_delay,1)];
data_h1 = [zeros(length(ramp_h2)+motion_delay,1); hprof1; zeros(length(ramp_h2)-motion_delay,1)];

% trailig motion profiles
data_p2 = [ramp_p2; pprof2+offset_p2; flip(ramp_p2)];
data_h2 = [ramp_h2; hprof2+offset_h2; flip(ramp_h2)];

% sinchronization command

% ref_sig = 
% data_sync = [zeros(length(ramp_p2)+motion_delay,1); ref_sig; zeros(length(ramp_p2)-motion_delay,1)];
% data_sync(end) = 0;

% time vector to create time series objects
toime = (0:length(data_p2)-1)'/srate;

% creating time series to be used by Simulink model
outp1 = timeseries(data_p1,toime);
outh1 = timeseries(data_h1,toime);
outp2 = timeseries(data_p2,toime);
outh2 = timeseries(data_h2,toime);
sync_sig = timeseries(data_sync,toime);

%% Plot the resulting profiles

figure()
subplot(3,1,1)
plot(toime,data_p1); hold on;
plot(toime,data_p2); hold off;
legend('pitch1', 'pitch2')
subplot(3,1,2)
plot(toime,data_h1); hold on;
plot(toime,data_h2); hold off;
legend('heave1', 'heave2')
subplot(3,1,3)
plot(toime,data_sync);
ylim([0,2])
legend('reference')

%% For zero motion

% outp1 = timeseries(0,1);
% outh1 = timeseries(0,1);
% 
% outp2 = timeseries(0,1);
% outh2 = timeseries(0,1);

%% Required simulation time

sim_time = ceil(toime(end))+5

%% Run here

running = 1; % matlab code pause variable
% set_param('simulink_traverse_control','SimulationCommand','start');
set_param('test_model','SimulationCommand','start');

disp('running traverse...')
pause(sim_time+10);

%% Convert output data from Simulink

F_wallace = conv_output_wallace_force(force_wallace,0);
F_gromit = conv_output_gromit_force(force_gromit,0);

encoders = [encoder_pitch1, encoder_heave1, encoder_pitch2, encoder_heave2];
encoders(encoders(:,1:4)>1e6) = encoders(encoders(:,1:4)>1e6)-2^32;
encoders(:,1) = (encoders(:,1)).*2*pi/10000;

