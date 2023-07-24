%% Load data, see Libraries/Analysis/dataLocations.m for more data storage location strings
datadir = 'D:\Experiments\2foil\SmallFoilAndVib_pitch1=20deg,heave1=0.5c';
thcknss = 0.0265; span_foil =0.365;%chord_foil = 0.06; %cross
%% Find the data files
dir_fullname = [datadir,'\data\'];
trialfiles = dir([dir_fullname,'bias*.mat']);
% trialfiles=trialfiles(~contains({trialfiles.name},'bias')); 
trialfiles=natsortfiles(trialfiles);
numtrials = length(trialfiles);

%% initialize arrays
trialstep=1;
biasWallace = nan(floor(numtrials/trialstep),6);

%% Load bias files
i=1;
for trial_index = 1:trialstep:numtrials
    try
        load([dir_fullname,trialfiles(trial_index).name]);
    catch
        disp(['Failed to load ',trialname])
    end
    biasWallace(i,:) = bias.Wallace;
    i=i+1;
end

plot(biasWallace./mean(biasWallace))