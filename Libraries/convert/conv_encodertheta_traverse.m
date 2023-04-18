function positionDeg = conv_encodertheta_traverse(encoderData, positionOffset)
%[positionDeg] = conv_encodertheta_traverse(encoderData)
%
% Belt traverse theta/pitch axis encoder data convertor.
% Converts counter channel measurements [counts] to relative rotation [deg].
%
% v1.0.0, xh, 02/01/2023
% -------------------------------------------------------------------------
% Inputs:
% encoderData = encoder measurement from the DAQ counter in [counts],
%               one-dimensional vector;
%
% positionOffset = initial/reference position of the measurement in [deg].
% -------------------------------------------------------------------------
% Output:
% positionDeg = position output referenced to postionOffset in [deg].
% -------------------------------------------------------------------------
% Note:
% Change encoderCPR, encoderType, and counterNBits accordingly if
% parameters of the encoder and the bit-rate of the DAQ change.
% -------------------------------------------------------------------------

encoderCPR = 2000; % counts per revolution = 2000
encoderType = 4; % encoder type according to the DAQ counter channel setup

counterNBits = 32; % counter channel resolution [32-bit]
signedThreshold = 2^(counterNBits-1); % count range
% obtain signed count values
encoderData(encoderData>signedThreshold) = encoderData(encoderData>signedThreshold) - 2^counterNBits;

% convert counts to [deg] relative to the initial offset position
% theta axis pointing downward, positive is clockwise looking from above
positionDeg = - encoderData*360/(encoderType*encoderCPR) + positionOffset;

end