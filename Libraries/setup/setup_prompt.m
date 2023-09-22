function [Parameters] = setup_prompt()
% default answers
P.sampleRate = '1000';
P.firstFoilShape = 'V1';
P.secondFoilShape = 'None';
P.flumeDepthMeters = '0.55';
P.flumeHertz = '16.0';
P.foilSeparationMeters = '0.242'; 
P.pitchOffsetDegG = '0';
P.heaveOffsetMetersG = '0';
P.pitchOffsetDegW = '181';
P.heaveOffsetMetersW = '0.23';
P.temperatureCelsius = '22.77';
P.pitchAxis = '0.5';
P.pivFlag = 'n';
P.filterFlag = 'n';
P.experimentName = 'Enter descriptive name';
P.save2LRS = 'y';
P.wallDistanceLeftMeters = '0.4';
P.wallDistanceRightMeters = '0.4';
% expected time delay between Gromit and Wallace (Gromit leading motion)
P.pitchDelayW = 13;
P.heaveDelayW = 50; % phase difference between pitch wallace and heave wallace (heave lags behind pitch)
disp(['NOTE: Expected time delay between Gromit and Wallace motions (Gromit leading the motion) is set to ', ...
    num2str(P.pitchDelayW),' ms']);

% generate prompt window
prompt = {'Enter the sample rate (in Hertz)','Enter first foil shape (as string): ','Enter second foil shape (as string): ', ...
    'Enter Flume water depth (in meters): ', 'Enter anticipated flume frequency (Hz): ',  'Enter foil separation distance (m): ', ...
    'Enter Gromit pitch offset from starting position (deg): ','Enter Gromit heave offset from starting position (m): ',...
    'Enter Wallace pitch offset from homing position (deg): ','Enter Wallace heave offset from homing position (m): ',...
    'Enter flume water temperature (from vectrino, in degC):   ','Enter foil pitch axis position (0 to 1 chord)', ...
    'Using PIV? (y/n)', 'Filter data at 60Hz? (y/n)','Enter Experiment Folder Name:','Save to LRS? (y/n):',...
    'Enter Mean wall distance (left, in meters): ','Enter Mean wall distance (right, in meters): '};

% input dialog window
numLines = 1; 
expCellArray = struct2cell(P);
answer = inputdlg(prompt, 'Experiment configuration', numLines, expCellArray);

% store reponses
P.sampleRate = str2double(answer{1});
P.firstFoilShape = char(answer{2});
P.secondFoilShape = char(answer{3});
P.flumeDepthMeters = str2double(answer{4});
P.flumeHertz = str2double(answer{5});
P.foilSeparationMeters = str2double(answer{6});
P.pitchOffsetDegG = str2double(answer{7});
P.heaveOffsetMetersG = str2double(answer{8});
P.pitchOffsetDegW = str2double(answer{9});
P.heaveOffsetMetersW = str2double(answer{10});
P.temperatureCelsius = str2double(answer{11});
P.pitchAxis = str2double(answer{12});
P.pivFlag = char(answer{13});
P.filterFlag = char(answer{14});
P.experimentName = char(answer{15});
P.save2LRS = char(answer{16});
P.wallDistanceLeftMeters = str2double(answer{17});
P.wallDistanceRightMeters = str2double(answer{18});

 % finds properties of the foils selected for the experiment
P.Foils = foils_database(P.firstFoilShape,P.secondFoilShape);

% establish save folder
if P.save2LRS == 'y'
    P.dataFolderName = ['R:\ENG_Breuer_Shared\group\JoelNewbolt\ExperimentalData\',P.experimentName];
else
    P.dataFolderName = ['D:\Experiments\Data\',P.experimentName];
end

% Prompt user for a new folder name if the previous folder already exists (overwrite protection)
[~, msg, ~] = mkdir(P.dataFolderName);
while strcmp(msg,'Directory already exists.')
    disp('Tried to create data folder that already exitsts! Appending date and time to folder name.')
    currDate = strrep(strrep(char(datetime), ':', '-'),' ','_');
    newDataDir = [P.experimentName,'_',currDate];
    if P.save2LRS == 'y'
        P.dataFolderName = ['R:\ENG_Breuer_Shared\group\JoelNewbolt\ExperimentalData\',newDataDir];
    else
        P.dataFolderName = ['D:\Experiments\Data\',newDataDir];
    end
    [~, msg, ~] = mkdir(P.dataFolderName);
end
folder_name = [P.dataFolderName,'\data'];
mkdir(folder_name); % for the bias measurements

% e.offset_home = [e.firstFoilPitchOffsetDegrees, e.firstFoilHeaveOffsetMeters, e.secondFoilPitchOffsetDegrees, e.secondFoilHeaveOffsetMeters];
Parameters = P;
end
