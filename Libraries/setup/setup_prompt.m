function [ExperimentParameters] = setup_prompt()
% default answers
EP.sampleRate = '1000';
EP.firstFoilShape = 'V1';
EP.secondFoilShape = 'None';
EP.flumeDepthMeters = '0.55';
EP.flumeHertz = '16.0';
EP.foilSeparationMeters = '0.242'; 
EP.pitchOffsetDegG = '0';
EP.heaveOffsetMetersG = '0';
EP.pitchOffsetDegW = '181';
EP.heaveOffsetMetersW = '0.23';
EP.temperatureCelsius = '22.77';
EP.pitchAxis = '0.5';
EP.pivFlag = 'n';
EP.filterFlag = 'n';
EP.experimentName = 'Enter descriptive name';
EP.save2LRS = 'y';
EP.wallDistanceLeftMeters = '0.4';
EP.wallDistanceRightMeters = '0.4';
% expected time delay between Gromit and Wallace (Gromit leading motion)
EP.pitchDelayW = 13;
EP.heaveDelayW = 50; % phase difference between pitch wallace and heave wallace (heave lags behind pitch)
disp(['NOTE: Expected time delay between Gromit and Wallace motions (Gromit leading the motion) is set to ', ...
    num2str(EP.pitchDelayW),' ms']);

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
expCellArray = struct2cell(EP);
answer = inputdlg(prompt, 'Experiment configuration', numLines, expCellArray);

% store reponses
EP.sampleRate = str2double(answer{1});
EP.firstFoilShape = char(answer{2});
EP.secondFoilShape = char(answer{3});
EP.flumeDepthMeters = str2double(answer{4});
EP.flumeHertz = str2double(answer{5});
EP.foilSeparationMeters = str2double(answer{6});
EP.pitchOffsetDegG = str2double(answer{7});
EP.heaveOffsetMetersG = str2double(answer{8});
EP.pitchOffsetDegW = str2double(answer{9});
EP.heaveOffsetMetersW = str2double(answer{10});
EP.temperatureCelsius = str2double(answer{11});
EP.pitchAxis = str2double(answer{12});
EP.pivFlag = char(answer{13});
EP.filterFlag = char(answer{14});
EP.experimentName = char(answer{15});
EP.save2LRS = char(answer{16});
EP.wallDistanceLeftMeters = str2double(answer{17});
EP.wallDistanceRightMeters = str2double(answer{18});

 % finds properties of the foils selected for the experiment
EP.Foils = foils_database(EP.firstFoilShape,EP.secondFoilShape);

% establish save folder
if EP.save2LRS == 'y'
    EP.dataFolderName = ['R:\ENG_Breuer_Shared\group\JoelNewbolt\ExperimentalData\',EP.experimentName];
else
    EP.dataFolderName = ['D:\Experiments\Data\',EP.experimentName];
end

% Prompt user for a new folder name if the previous folder already exists (overwrite protection)
[~, msg, ~] = mkdir(EP.dataFolderName);
while strcmp(msg,'Directory already exists.')
    disp('Tried to create data folder that already exitsts! Appending date and time to folder name.')
    currDate = strrep(strrep(char(datetime), ':', '-'),' ','_');
    newDataDir = [EP.experimentName,'_',currDate];
    if EP.save2LRS == 'y'
        EP.dataFolderName = ['R:\ENG_Breuer_Shared\group\JoelNewbolt\ExperimentalData\',newDataDir];
    else
        EP.dataFolderName = ['D:\Experiments\Data\',newDataDir];
    end
    [~, msg, ~] = mkdir(EP.dataFolderName);
end
folder_name = [EP.dataFolderName,'\data'];
mkdir(folder_name); % for the bias measurements

% e.offset_home = [e.firstFoilPitchOffsetDegrees, e.firstFoilHeaveOffsetMeters, e.secondFoilPitchOffsetDegrees, e.secondFoilHeaveOffsetMeters];
ExperimentParameters = EP;
end
