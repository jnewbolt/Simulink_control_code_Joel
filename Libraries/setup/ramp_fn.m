function [timeRamp, ramp] = ramp_fn(rampTime, startPosition, endPosition, sampleRate)
    % used by both traverses in order to ramp up and down from starting position to offset and bias
    numPtsRamp = rampTime*sampleRate;
    timeRamp = (0:numPtsRamp-1)/sampleRate; % dimensional ramping time vector
    
    ramp = startPosition + (endPosition-startPosition)*(0.5*(1-cos( pi*(0:numPtsRamp-1)/numPtsRamp)))';
end