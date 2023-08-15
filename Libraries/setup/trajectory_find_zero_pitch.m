function [times, pitchDegreesG, heaveMetersG, pitchDegreesW, heaveMetersW, syncSig] = trajectory_find_zero_pitch(EP, scanTime, ...
    pitchAmpDeg, traverse)
        
    T = 1/EP.sampleRate; % define timestep in seconds

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
            profs(:,1) = pprof1+EP.firstFoilPitchOffsetDegrees;
            profs(:,2) = ones(size(pprof1))*EP.firstFoilHeaveOffsetMeters;
            profs(:,3) = zeros(size(pprof1))+EP.secondFoilPitchOffsetDegrees;
        case 'Wallace'
            profs(:,3) = pprof1+EP.secondFoilPitchOffsetDegrees;
            profs(:,1) = zeros(size(pprof1))+EP.firstFoilPitchOffsetDegrees;
            % Move Gromit out of the way to heave position -0.15 meters
            numPtsRamp = length(pprof1_a);
            offsetMetersG = (EP.firstFoilHeaveOffsetMeters-0.15);
            rampHeaveMetersG = offsetMetersG*(0.5*(1-cos( pi*(0:numPtsRamp-1)/numPtsRamp)))';
            profs(:,2) = [rampHeaveMetersG; offsetMetersG*ones(length(pprof1)-2*length(pprof1_a),1); flip(rampHeaveMetersG)];
    end
    profs(:,4) = ones(size(pprof1))*EP.secondFoilHeaveOffsetMeters;
    profs(:,5) = zeros(size(pprof1));
    % mark the locations at which data will be extracted for alignment calculations:
    profs([stp1, stp1+stp2, stp1+stp2+stp3, stp1+stp2+stp3+stp4],5) = 1;

    
    % convert into time series to be output to simulink
    times = (0:size(profs,1)-1)'/EP.sampleRate; % time vector to create time series objects
    pitchDegreesG = timeseries(profs(:,1),times);
    heaveMetersG = timeseries(profs(:,2),times);
    pitchDegreesW = timeseries(profs(:,3),times);
    heaveMetersW = timeseries(profs(:,4),times);
    syncSig = timeseries(profs(:,5),times);

    % plot trajectories
    plot_profiles(pitchDegreesG,heaveMetersG,pitchDegreesW,heaveMetersW,syncSig);

end