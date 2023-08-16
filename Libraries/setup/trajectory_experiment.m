function [Time, Waveform] = trajectory_experiment(nCycles, frequency, sampleRate, nCycsRampUp, nCycsRampDown, amplitude, phaseLag, constantAmplitude)
% Function to generate a waveform of a given amplitude, with a given frequency [Hz], sampling
% frequency [Hz], number of cycles to ramp up and down, Phase shift [deg]

%%
dt = 1/sampleRate; % time steps for the sine wave
Npts_total = round((nCycles + nCycsRampUp + nCycsRampDown)*(1/frequency)*sampleRate); % total number of points
                                                                                      % (cycles + transients) x period x sampling frequency
Npts_per_cycle = round(Npts_total/(nCycles+nCycsRampUp+nCycsRampDown)); % total number of points per cycle

Time = (0:Npts_total)*dt; % time vector

% Base waveform - the begining and end will be modulated
Waveform = sin(2*pi*frequency*Time - deg2rad(phaseLag));
if constantAmplitude ~= 0
    Waveform = ones(1, length(Time));
end
%% Now add in the ramp up and ramp down

% Modulate the start and the end with a cosine window:

Npts_rampup = round(Npts_per_cycle*nCycsRampUp); % number of cycles to be used in the transient
Ramp_up_Function = 0.5*(1-cos( pi*(0:Npts_rampup-1)/Npts_rampup));
index_up = 1:Npts_rampup;
Waveform(index_up) = Ramp_up_Function .* Waveform(index_up); % replaces cycles of transient with the base wave x rampup function

Npts_rampdown = round(Npts_per_cycle*nCycsRampDown); % same process but for the last cycles
Ramp_down_Function = 0.5*(1+cos( pi*(0:Npts_rampdown-1)/Npts_rampdown));
% index_down = round(Npts_total-Npts_rampdown+1):Npts_total; % needs to be a rounded value in order to be indexable
index_down = Npts_total-Npts_rampdown+1:Npts_total;
Waveform(index_down) = Ramp_down_Function .* Waveform(index_down);

Waveform(end) = 0; % probably when rounding the indeces, it rounds down the last index and therefore leaves it out from the ramp-down,
                   % therefore we can choose the last point to be zero

%% Multiply by the amplitude

Waveform = (amplitude*Waveform)';

end
