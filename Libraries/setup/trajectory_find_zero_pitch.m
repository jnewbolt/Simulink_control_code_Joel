function [times, pitchDegreesG, heaveMetersG, pitchDegreesW, heaveMetersW, syncSig] = trajectory_find_zero_pitch(EP, scanTime, ...
    pitchAmpDeg, traverse)
        
    T = 1/EP.sampleRate; % define timestep in seconds
    switch traverse
        case 'Wallace'
            EP.firstFoilHeaveOffsetMeters = -0.15; % Gromit moves out of the way
    end
    % ramp from home position for gromit traverse
    rampTime = 5; % Time to move to starting position, in seconds
    [~, ramp_p1, ramp_h1] = ramp_fn(rampTime, EP, 'Gromit');
    [~, ramp_p2, ramp_h2] = ramp_fn(rampTime, EP, 'Wallace');
    
    % pitch motion profile will be constructed in steps:
    % a: move from zero to starting position
    pprof1_a = linspace(0, pitchAmpDeg, 0.5*scanTime/T)';
    % b: wait at starting position (stp1 is to be used by the reference signal)
    pprof1_b = linspace(pitchAmpDeg, pitchAmpDeg, 0.5*scanTime/T)'; stp1 = length(pprof1_a)+length(pprof1_b);
    % c: move to opposite position (stp2 is to be used by the reference signal)
    pprof1_c = linspace(pitchAmpDeg, -pitchAmpDeg, scanTime/T)'; stp2 = length(pprof1_c);
    % d: wait at opposite position (stp3 is to be used by the reference signal)
    pprof1_d = linspace(-pitchAmpDeg, -pitchAmpDeg, scanTime/T)'; stp3 = length(pprof1_d);
    % e: move back to starting position (stp4 is to be used by the reference signal)
    pprof1_e = linspace(-pitchAmpDeg, pitchAmpDeg, scanTime/T)'; stp4 = length(pprof1_e);
    % f: wait at starting position
    pprof1_f = linspace(pitchAmpDeg, pitchAmpDeg, 0.5*scanTime/T)';
    % g: move to zero
    pprof1_g = linspace(pitchAmpDeg, 0, 0.5*scanTime/T)';
    
    % concatenate profiles into a single commanded profile
    pprof1 = [pprof1_a; pprof1_b; pprof1_c; pprof1_d; pprof1_e; pprof1_f; pprof1_g];
    
    % concatenate ramps with commanded motion profiles
    switch traverse
        case 'Gromit'
            profs(:,1) = [ramp_p1; pprof1+EP.firstFoilPitchOffsetDegrees; flip(ramp_p1)];
            profs(:,3) = zeros(size(profs(:,1)));
        case 'Wallace'
            profs(:,3) = [ramp_p2; pprof1+EP.secondFoilPitchOffsetDegrees; flip(ramp_p2)];
            profs(:,1) = zeros(size(profs(:,1)));
    end
    profs(:,2) = [ramp_h1; ones(size(pprof1))*EP.firstFoilHeaveOffsetMeters; flip(ramp_h1)];
    profs(:,4) = [ramp_h2; ones(size(pprof1))*EP.secondFoilHeaveOffsetMeters; flip(ramp_h2)];
    profs(:,5) = [zeros(size(ramp_p1)); zeros(size(pprof1)); zeros(size(ramp_p1))];
    % mark the locations at which data will be extracted for alignment calculations:
    profs([stp1, stp1+stp2, stp1+stp2+stp3, stp1+stp2+stp3+stp4]+length(ramp_p1),5) = 1;

    
    % convert into time series to be output to simulink
    times = (0:size(profs,1)-1)'/EP.sampleRate; % time vector to create time series objects
    pitchDegreesG = timeseries(profs(:,1),times);
    heaveMetersG = timeseries(profs(:,2),times);
    pitchDegreesW = timeseries(profs(:,3),times);
    heaveMetersW = timeseries(profs(:,4),times);
    syncSig = timeseries(profs(:,5),times);

    switch traverse
        case 'Gromit'
            pitchDegreesW = timeseries(zeros(size(times)),times); % Zero motion trajectory for wallace
            heaveMetersW = timeseries(zeros(size(times)),times);
    end

    % plot trajectories
    plot_profiles(pitchDegreesG,heaveMetersG,pitchDegreesW,heaveMetersW,syncSig);

end