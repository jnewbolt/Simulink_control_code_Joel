%% Find zero pitch simulink

% function bias = find_zero_pitch_simulink(experiment, bias, traverse)

stndby = 20;

offset_home = experiment.offset_home;

switch traverse
    case 'g' % leading traverse (Gromit)
        fy_index = 18; % index of the perpendicular (to the streamwise direction) transducer force channel
        mz_index = 22; % index of the pitch axis transducer moment channel
        prof_index = 1; % index of the motion profile corresponding to the selected traverse
        p_bias_index = 1; % pitch bias index for the corresponding traverse
        traverse_name = 'Gromit';
    case 'w' % trailing traverse (Wallace)
        fy_index = 8;
        mz_index = 12;
        prof_index = 3;
        p_bias_index = 2;
        traverse_name = 'Wallace';
    otherwise
        error('Unexpected traverse name. Argument invalid.')
end

%% Alignment

disp(['Aligning ',traverse_name, ' traverse.']);

repeat_alignment = 1;
while repeat_alignment == 1
    %% Generate profiles
    % for initial search:
    scan_time = 5;
    search_amplitude = 5;
    expand_search = 0;

    while expand_search < 2
        expand_search = expand_search + 1;

        [toime, outp1, outh1, outp2, outh2, sync_sig] = alingment_profile(experiment, traverse, scan_time, search_amplitude, bias);
    
        % simulation time
        sim_time = ceil(toime(end))+2;
        disp(['Expected simulation time: ', num2str(sim_time), ' seconds']);
        
        % pass parameters for gromit heave gain in simulation
        freq = 0; heave1 = 0;
    
        %% Run traverse

        % clear variables before next experiment
        clear raw_encoder_p1 raw_encoder_h1 raw_encoder_p2 raw_encoder_h2 raw_force_wallace raw_force_gromit ref_signal
    
        disp(['Finding zero.  Traverse will now move to +/- ',num2str(search_amplitude),' degrees']);
        set_param('simulink_traverse_control','SimulationCommand','start');
        
        disp('Running traverse...')
        % NOTE: this needs to account for the build time
        pause(sim_time);
        disp('Acquiring data...')
        
        % check if data is
        while ~exist('raw_encoder_p1','var') || ~exist('raw_encoder_h1','var') || ~exist('raw_encoder_p2','var') || ~exist('raw_encoder_h2','var') || ~exist('raw_force_wallace','var') || ~exist('raw_force_gromit','var') || ~exist('ref_signal','var')
            pause(5)
            disp('Loading...')
        end

        disp('Done')
    
        %% Convert data
    
        raw_encoders = [raw_encoder_p1, raw_encoder_h1, raw_encoder_p2, raw_encoder_h2];
        out = convert_output(raw_encoders, raw_force_wallace, raw_force_gromit, raw_vectrino, ref_signal, bias, [], offset_home);
    
        %% Extract relevant forces for analysis
    
        % inices of non-zero elements
        k = find(ref_signal == 1);
        
        range_pos = k(1):k(2); % range of positive angle sweep
        range_neg = k(3):k(4); % range of negative angle sweep
    
        Fy_pos = out(range_pos,fy_index);
        Ty_pos = out(range_pos,mz_index);
        pitch_pos = out(range_pos,prof_index);
    
        Fy_neg = out(range_neg,fy_index);
        Ty_neg = out(range_neg,mz_index);
        pitch_neg = out(range_neg,prof_index);
    
        b_Vtheta_pos = pitch_pos(1);
        a_Vtheta_pos = (pitch_pos(end)-b_Vtheta_pos)/(scan_time*experiment.srate);
    
        b_Vtheta_neg = pitch_neg(1);
        a_Vtheta_neg = (pitch_neg(end)-b_Vtheta_neg)/(scan_time*experiment.srate);
    
        %% Linear fits on curves
    
        coefL1_pos = polyfit(1:numel(Fy_pos), smooth(Fy_pos,100)', 1);
        aL_pos = coefL1_pos(1);
        bL_pos = coefL1_pos(2);
    
        coefL1_neg = polyfit(1:numel(Fy_neg), smooth(Fy_neg,100)', 1);
        aL_neg = coefL1_neg(1);
        bL_neg = coefL1_neg(2);
    
        plot(1:numel(pitch_pos), pitch_pos, 1:numel(Fy_pos), smooth(Fy_pos,100), 1:numel(Fy_pos), aL_pos*(1:numel(Fy_pos))+bL_pos); hold on;
        plot(1:numel(pitch_neg), pitch_neg, 1:numel(Fy_neg), smooth(Fy_neg,100), 1:numel(Fy_neg), aL_neg*(1:numel(Fy_neg))+bL_neg); hold off;
    
        legend('Pitch negative','Fy smoothed','Linear fit','Pitch positive','Fy smoothed','Linear fit');
    
        %% Expanding search
    
        if max(aL_pos*(1:numel(Fy_pos))+bL_pos) < 0 || min(aL_pos*(1:numel(Fy_pos))+bL_pos) > 0
            if expand_search == 2
                disp('Pitch out of range. Align manually.')
                break
            else
                disp('Pitch out of range. Expanding search...')
                scan_time = 8;
                search_amplitude = 10;
            end
        else
            break
        end

    end % expand search end

    %% Calculate the average pitch bias

    pitch_bias_pos = a_Vtheta_pos*(-bL_pos/aL_pos) + b_Vtheta_pos;
    pitch_bias_neg = a_Vtheta_neg*(-bL_neg/aL_neg) + b_Vtheta_neg;
    new_bias = mean([pitch_bias_pos, pitch_bias_neg]);

    disp(['Pitch Bias [deg]: ', num2str(new_bias)])

%     pitch_check = input(['Does this look alright to you? y/n + Enter', newline],"s");
    % bypass input
    pitch_check = 'y';

    if strcmp(pitch_check,'y')
        bias.pitch(p_bias_index) = new_bias;
        disp([traverse_name, ' pitch bias updated'])
        break; % terminate while loop
    else
        disp([traverse_name, ' pitch bias not updated'])
        repeat_response = input(['Would you like to repeat find_zero_pitch for ',traverse_name,'? y/n + Enter',newline],"s");
        if strcmp(repeat_response,'n')
            break; % terminate while loop
        end
    end

end % repeat alignment search

disp([traverse_name, ' alignment done.']);


%% In-routine functions

    function [toime, outp1, outh1, outp2, outh2, sync_sig] = alingment_profile(experiment, traverse, t, ppos, bias)
        
        T = 1/experiment.srate;
        
        switch traverse

            case 'g'
                offset = experiment.offset_home;
                % no additional offset during alignment
                offset(3) = 0; % wallace stays at homing position
                offset(4) = 0; % wallace stays at homing position
                % ramp from home position for gromit traverse
                [~, ramp_p1, ramp_h1] = ramp_fn(5, T, bias, offset, 'g');
                [~, ramp_p2, ramp_h2] = ramp_fn(5, T, bias, offset, 'w');

                % pitch motion profile will be constructed in steps:
                % a: move from zero to starting position
                pprof1_a = linspace(0, ppos, 0.5*t/T)';
                % b: wait at starting position (stp1 is to be used by the reference signal)
                pprof1_b = linspace(ppos, ppos, 0.5*t/T)'; stp1 = length(pprof1_a)+length(pprof1_b);
                % c: move to opposite position (stp2 is to be used by the reference signal)
                pprof1_c = linspace(ppos, -ppos, t/T)'; stp2 = length(pprof1_c);
                % d: wait at opposite position (stp3 is to be used by the reference signal)
                pprof1_d = linspace(-ppos, -ppos, t/T)'; stp3 = length(pprof1_d);
                % e: move back to starting position (stp4 is to be used by the reference signal)
                pprof1_e = linspace(-ppos, ppos, t/T)'; stp4 = length(pprof1_e);
                % f: wait at starting position
                pprof1_f = linspace(ppos, ppos, 0.5*t/T)';
                % g: move to zero
                pprof1_g = linspace(ppos, 0, 0.5*t/T)';
                
                % concatenate profiles into a single commanded profile
                pprof1 = [pprof1_a; pprof1_b; pprof1_c; pprof1_d; pprof1_e; pprof1_f; pprof1_g];
                % heave will not move for the whole trajectory of pitch
                hprof1 = zeros(size(pprof1));
                
                % zero motion trajectory for wallace
                pprof2 = zeros(size(pprof1));
                hprof2 = zeros(size(pprof1));

                % concatenate ramps with commanded motion profiles
                profs(:,1) = [ramp_p1; pprof1+offset(1); flip(ramp_p1)];
                profs(:,2) = [ramp_h1; hprof1+offset(2); flip(ramp_h1)];
                profs(:,3) = [ramp_p2; pprof2+offset(3); flip(ramp_p2)];
                profs(:,4) = [ramp_h2; hprof2+offset(4); flip(ramp_h2)];
                profs(:,5) = [zeros(size(ramp_p1)); zeros(size(pprof1)); zeros(size(ramp_p1))];
                % mark the locations at which data will be extracted for alignment calculations:
                profs([stp1, stp1+stp2, stp1+stp2+stp3, stp1+stp2+stp3+stp4]+length(ramp_p1),5) = 1;

                % plot trajectories
                plot_profiles(profs);
                
                % convert into time series to be output to simulink
                toime = (0:size(profs,1)-1)'/experiment.srate; % time vector to create time series objects
                outp1 = timeseries(profs(:,1),toime);
                outh1 = timeseries(profs(:,2),toime);
                outp2 = timeseries(profs(:,3),toime);
                outh2 = timeseries(profs(:,4),toime);
                sync_sig = timeseries(profs(:,5),toime);

            case 'w'
                offset = experiment.offset_home;
                % no additional offset during alignment
                offset(1) = 0;
                offset(2) = -0.15; % gromit moves out of the way
                % ramp from home position for gromit traverse
                [~, ramp_p1, ramp_h1] = ramp_fn(5, T, bias, offset, 'g');
                [~, ramp_p2, ramp_h2] = ramp_fn(5, T, bias, offset, 'w');
                
                % pitch motion profile will be constructed in steps:
                % a: move from zero to starting position
                pprof2_a = linspace(0, ppos, 0.5*t/T)';
                % b: wait at starting position (stp1 is to be used by the reference signal)
                pprof2_b = linspace(ppos, ppos, 0.5*t/T)'; stp1 = length(pprof2_a)+length(pprof2_b);
                % c: move to opposite position (stp2 is to be used by the reference signal)
                pprof2_c = linspace(ppos, -ppos, t/T)'; stp2 = length(pprof2_c);
                % d: wait at opposite position (stp3 is to be used by the reference signal)
                pprof2_d = linspace(-ppos, -ppos, t/T)'; stp3 = length(pprof2_d);
                % e: move back to starting position (stp4 is to be used by the reference signal)
                pprof2_e = linspace(-ppos, ppos, t/T)'; stp4 = length(pprof2_e);
                % d: wait at starting position
                pprof2_f = linspace(ppos, ppos, 0.5*t/T)';
                % f: move to zero
                pprof2_g = linspace(ppos, 0, 0.5*t/T)';
                
                % concatenate profiles into a single commanded profile
                pprof2 = [pprof2_a; pprof2_b; pprof2_c; pprof2_d; pprof2_e; pprof2_f; pprof2_g];
                % heave will not move for the whole trajectory of pitch
                hprof2 = zeros(size(pprof2));
                
                % zero motion trajectory for gromit
                pprof1 = zeros(size(pprof2));
                hprof1 = zeros(size(pprof2));

                % concatenate ramps with commanded motion profiles
                profs(:,1) = [ramp_p1; pprof1+offset(1); flip(ramp_p1)];
                profs(:,2) = [ramp_h1; hprof1+offset(2); flip(ramp_h1)];
                profs(:,3) = [ramp_p2; pprof2+offset(3); flip(ramp_p2)];
                profs(:,4) = [ramp_h2; hprof2+offset(4); flip(ramp_h2)];
                profs(:,5) = [zeros(size(ramp_p1)); zeros(size(pprof1)); zeros(size(ramp_p1))];
                % mark the locations at which data will be extracted for alignment calculations:
                profs([stp1, stp1+stp2, stp1+stp2+stp3, stp1+stp2+stp3+stp4]+length(ramp_p1),5) = 1;

                % plot trajectories
                plot_profiles(profs);
                
                % convert into time series to be output to simulink
                toime = (0:size(profs,1)-1)'/experiment.srate; % time vector to create time series objects
                outp1 = timeseries(profs(:,1),toime);
                outh1 = timeseries(profs(:,2),toime);
                outp2 = timeseries(profs(:,3),toime);
                outh2 = timeseries(profs(:,4),toime);
                sync_sig = timeseries(profs(:,5),toime);
        end
    end


% end