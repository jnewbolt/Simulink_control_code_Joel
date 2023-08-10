function [t_ramp, ramp_p, ramp_h] = ramp_fn(rampTime, ExperimentParameters, traverse)
    % used by both traverses in order to ramp up and down from starting position to offset and bias
    EP = ExperimentParameters;
    switch traverse
        case 'g'
            offsetPitch = EP.firstFoilPitchOffsetDegrees;
            offsetHeave = EP.firstFoilHeaveOffsetMeters;
        case 'w'
            offsetPitch = EP.secondFoilPitchOffsetDegrees;
            offsetHeave = EP.secondFoilHeaveOffsetMeters;
    end

    numPtsRamp = rampTime*EP.sampleRate;
    t_ramp = (0:numPtsRamp-1)/EP.sampleRate; % dimensional ramping time vector
    
    ramp_p = (offsetPitch)*(0.5*(1-cos( pi*(0:numPtsRamp-1)/numPtsRamp)))';
    ramp_h = (offsetHeave)*(0.5*(1-cos( pi*(0:numPtsRamp-1)/numPtsRamp)))';

end