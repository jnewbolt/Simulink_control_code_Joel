function positionm = conv_encodery_traverse(encoderData, positionOffset)
%[positionm] = conv_encodery_traverse(encoderData)
%
% Belt traverse y/heave axis encoder data convertor.
% Converts counter channel measurements [counts] to relative rotation [deg].
%
% v1.0.0, xh, 02/01/2023
% -------------------------------------------------------------------------
% Inputs:
% encoderData = encoder measurement from the DAQ counter in [counts],
%               one-dimensional vector;
%
% positionOffset = initial/reference position of the measurement in [m].
% -------------------------------------------------------------------------
% Output:
% positionm = position output referenced to postionOffset in [m].
% -------------------------------------------------------------------------
% Note:
% Change encoderCPR, encoderType, and counterNBits accordingly if
% parameters of the encoder and the bit-rate of the DAQ change.
% -------------------------------------------------------------------------

encoderCPR = 40; % counts per 2 mm = 40
encoderType = 4; % encoder type according to the DAQ counter channel setup

counterNBits = 32; % counter channel resolution [32-bit]
signedThreshold = 2^(counterNBits-1); % count range
% obtain signed count values
encoderData(encoderData>signedThreshold) = encoderData(encoderData>signedThreshold) - 2^counterNBits;

% convert counts to [deg] relative to the initial offset position
% theta axis pointing downward, positive is clockwise looking from above
% positionm = - encoderData*0.002/(encoderType*encoderCPR) + positionOffset;
positionm = - encoderData*0.002/encoderCPR + positionOffset;

end