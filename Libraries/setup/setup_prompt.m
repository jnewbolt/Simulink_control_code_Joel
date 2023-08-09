function [experiment] = setup_prompt()
% default answers
e.sampleRate = '1000';
e.firstFoilShape = 'V1';
e.secondFoilShape = 'None';
e.flumeDepthMeters = '0.55';
e.flumeHertz = '16.0';
e.foilSeparationMeters = '0.242'; 
e.firstFoilPitchOffsetDegrees = '0';
e.firstFoilHeaveOffsetMeters = '0';
e.secondFoilPitchOffsetDegrees = '181';
e.secondFoilHeaveOffsetMeters = '0.23';
e.temperatureCelsius = '22.77';
e.pitchAxis = '0.5';
e.pivFlag = 'n';
e.filterFlag = 'n';
e.experimentName = 'Enter descriptive name';
e.save2LRS = 'y';
e.wallDistanceLeftMeters = '0.4';
e.wallDistanceRightMeters = '0.4';

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
expCellArray = struct2cell(e);
answer = inputdlg(prompt, 'Experiment configuration', numLines, expCellArray);

% store reponses
e.sampleRate = str2double(answer{1});
e.firstFoilShape = char(answer{2});
e.secondFoilShape = char(answer{3});
e.flumeDepthMeters = str2double(answer{4});
e.flumeHertz = str2double(answer{5});
e.foilSeparationMeters = str2double(answer{6});
e.firstFoilPitchOffsetDegrees = str2double(answer{7});
e.firstFoilHeaveOffsetMeters = str2double(answer{8});
e.secondFoilPitchOffsetDegrees = str2double(answer{9});
e.secondFoilHeaveOffsetMeters = str2double(answer{10});
e.temperatureCelsius = str2double(answer{11});
e.pitchAxis = str2double(answer{12});
e.pivFlag = char(answer{13});
e.filterFlag = char(answer{14});
e.experimentName = char(answer{15});
e.save2LRS = char(answer{16});
e.wallDistanceLeftMeters = str2double(answer{17});
e.wallDistanceRightMeters = str2double(answer{18});

% establish save folder
if e.save2LRS == 'y'
    e.dataFolderName = ['R:\ENG_Breuer_Shared\group\JoelNewbolt\ExperimentalData\',e.experimentName];
else
    e.dataFolderName = ['D:\Experiments\Data\',e.experimentName];
end

% Prompt user for a new folder name if the previous folder already exists (overwrite protection)
[status, msg, msgID] = mkdir(e.dataFolderName);
while strcmp(msg,'Directory already exists.')
    new_dir = input('Chosen directory name already exists!  Please type a new folder name below, then hit enter. \n',"s");
    e.dataFolderName = ['D:\Experiments\',num2str(e.Number_of_foils),'foil\',new_dir];
    [status, msg, msgID] = mkdir(e.dataFolderName);
end
folder_name = [e.dataFolderName,'\data'];
mkdir(folder_name); % for the bias measurements

% e.offset_home = [e.firstFoilPitchOffsetDegrees, e.firstFoilHeaveOffsetMeters, e.secondFoilPitchOffsetDegrees, e.secondFoilHeaveOffsetMeters];
experiment = e;
end
