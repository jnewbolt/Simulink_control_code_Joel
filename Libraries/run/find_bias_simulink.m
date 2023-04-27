%% Find bias - simulink control

% function [out, bias] = find_bias_simulink(experiment, bias, motion_delay, offset_home)
    
    % Find bias in the force measurements
    % NOTE: only overwrites the force biases unless bias is not specified

    % bias is an optional argument

    %% Input arguments

    offset_home = experiment.offset_home;
    motion_delay = experiment.motion_delay;

    %% Initialize heave and pitch biases

    if ~exist('bias','var') || isempty(bias)
        bias.pitch = [0,0];
        bias.heave = [0,0];
    end
    
    %% Stationary command profiles
    
    % ramp to offset time
    ramp_time = 5; % in [s]

    % Generate ramp profiles
    [~, ramp_p1, ramp_h1] = ramp_fn(ramp_time, experiment.T, bias, offset_home, 'g');
    [~, ramp_p2, ramp_h2] = ramp_fn(ramp_time, experiment.T, bias, offset_home, 'w');
    
    % Zero motion profile
    prof_time = 20; % length of time of bias measurement
    phprof12 = zeros(prof_time/experiment.T,1);
    
    clear profs % in case this scripts is run after a previous one
    % gromit pitch
    profs(:,1) = [zeros(motion_delay,1); ramp_p1; phprof12+offset_home(1)+bias.pitch(1); flip(ramp_p1)];
    % gromit heave
    profs(:,2) = [zeros(motion_delay,1); ramp_h1; phprof12+offset_home(2)+bias.heave(1); flip(ramp_h1)];
    % wallace pitch
    profs(:,3) = [ramp_p2; phprof12+offset_home(3)+bias.pitch(2); flip(ramp_p2); zeros(motion_delay,1)];
    % wallace heave
    profs(:,4) = [ramp_h2; phprof12+offset_home(4)+bias.heave(2); flip(ramp_h2); zeros(motion_delay,1)];
    
    % reference signal - has 0 for ramps and delays and 1 for usable data
    profs(:,5) = [zeros(motion_delay,1); zeros(size(ramp_p1)); ones(size(phprof12)); zeros(size(ramp_p1))];
    % NOTE: the profile for sync_sig should be value 1 some time AFTER the ramp is finished,
    % and should go back to 0 some time BEFORE the 0 profile ends.

    % plot trajectories
    plot_profiles(profs);
    
    % convert into time series to be output to simulink
    toime = (0:size(profs,1)-1)'/experiment.srate; % time vector to create time series objects
    outp1 = timeseries(profs(:,1),toime);
    outh1 = timeseries(profs(:,2),toime);
    outp2 = timeseries(profs(:,3),toime);
    outh2 = timeseries(profs(:,4),toime);
    sync_sig = timeseries(profs(:,5),toime);
    
    % simulation time
    sim_time = ceil(toime(end))+2;
    disp(['Expected simulation time: ', num2str(sim_time), ' seconds']);
    
    % pass parameters for gromit heave gain in simulation
    freq = 0; heave1 = 0;

    %% Run traverse

    % clear variables in case done after an experiment
    clear raw_encoder_p1 raw_encoder_h1 raw_encoder_p2 raw_encoder_h2 raw_force_wallace raw_force_gromit ref_signal

    set_param('simulink_traverse_control','SimulationCommand','start');
    
    disp('Running traverse...')
    pause(sim_time);
    disp('Acquiring data...')

    while ~exist('raw_encoder_p1','var') || ~exist('raw_encoder_h1','var') || ~exist('raw_encoder_p2','var') || ~exist('raw_encoder_h2','var') || ~exist('raw_force_wallace','var') || ~exist('raw_force_gromit','var') || ~exist('ref_signal','var')
            pause(5)
            disp('Loading...')
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
    
    range_x = find(ref_signal);
    range_x = range_x(1000:end);

    bias.Wallace = mean(raw_force_wallace(range_x,:),1);
    bias.WallaceStdev = std(raw_force_wallace(range_x,:),1);

    bias.Gromit = mean(raw_force_gromit(range_x,:),1);
    bias.GromitStdev = std(raw_force_gromit(range_x,:),1);

    bias.accmeter = mean(raw_accelmeter(range_x,:),1);
    bias.accmeterStdev = std(raw_accelmeter(range_x,:),1);

    %% Convert output

    raw_encoders = [raw_encoder_p1, raw_encoder_h1, raw_encoder_p2, raw_encoder_h2];
    out = convert_output(raw_encoders, raw_force_wallace, raw_force_gromit, raw_vectrino, raw_accelmeter, ref_signal, bias, range_x, offset_home);

    %% Plot results

    figure()
    force_resolution = [1/32 1/32 1/16];
    torque_resolution = 1/528;

    subplot(2,1,2)
    plot(out(20:end-20,7:9)./force_resolution,'.')
    hold on 
    plot(out(20:end-20,10:12)/torque_resolution,'.')
    hold off
    title('Wallace (last)')
    ylabel('Forces and Torques (normalized by resolution)')
    legend('Fy','Fx','Fz','Ty','Tx','Tz')

    subplot(2,1,1)
    plot(out(20:end-20,17:19)./force_resolution,'.')
    hold on 
    plot(out(20:end-20,20:22)/torque_resolution,'.')
    hold off
    title('Gromit (middle)')
    ylabel('Forces and Torques (normalized by resolution)')
    legend('Fy','Fx','Fz','Ty','Tx','Tz')
    
    %% Raise warnings for potentially faulty measurements

    if bias.accmeter < 1.5 || bias.accmeter > 1.8
        disp(bias.accmeter)
        disp('Warning: Accelerometer voltage outside expected range. Try power cycling.')
    end

    bias.RMSEW = sqrt(mean((raw_force_wallace - repmat(bias.Wallace,numel(raw_force_wallace(:,1)),1)).^2));
    if sum(bias.RMSEW>[.15 .15 .3 .1 .1 .1])>0 
            disp(bias.RMSEW);
        disp('Warning: Wallace error signal above normal. Check wiring/ grounding.')
    end
    bias.RMSEG = sqrt(mean((raw_force_gromit - repmat(bias.Gromit,numel(raw_force_gromit(:,1)),1)).^2));
    if sum(bias.RMSEG>[.15 .15 .3 .1 .1 .1])>0 
        disp('Warning: Gromit error signal above normal. Check wiring/ grounding.')
    end

    Percent_fullrange_error = bias.RMSEW./[660 660 1980 60 60 60]*100;
    Percent_fullrange_errorG = bias.RMSEG./[660 660 1980 60 60 60]*100;

    %% Save bias into a subfolder within the active saving folder
    
    time = clock;
    fname = experiment.fname; % extract from the experimental setup parameters

    folder_name = [fname,'\data'];
    
    numfiles = dir([folder_name,'\bias*']);
    jj = numel(numfiles)+1;
    filename = [folder_name,'\bias_',num2str(jj)];
    
    save(filename,'bias','time','out');

% end

