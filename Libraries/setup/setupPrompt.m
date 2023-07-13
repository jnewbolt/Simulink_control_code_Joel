function experiment = setupPrompt(fs)

% default answers

experiment.chord = 0.054;
experiment.thcknss = 0.054;
experiment.span = 0.401;
experiment.foil_shape = 'V1';
experiment.Wall_distance_left = 0.4;
experiment.Wall_distance_right = 0.4;
experiment.flume_height = 0.55;
experiment.flume_hertz = 16.0;
experiment.Number_of_foils = 2;
experiment.foil_separation = 0.242; 
experiment.foil_offset = 0;
experiment.offset_p1 = 0;
experiment.offset_h1 = 0;
experiment.offset_p2 = 181;
experiment.offset_h2 = 0.235;
experiment.Temperature = 22.77;
experiment.pitch_axis = 0.5;
experiment.piv_var = 0;
experiment.filt_var = 0;
experiment.expf_name = 'Enter descriptive name';
experiment.save_lrs = 'n';

defaultanswers = {num2str(experiment.chord),...
    num2str(experiment.span),...
    experiment.foil_shape,...
    num2str(experiment.Wall_distance_left),...
    num2str(experiment.Wall_distance_right),...
    num2str(experiment.flume_height),...
    num2str(experiment.flume_hertz),...
    num2str(experiment.Number_of_foils),...
    num2str(experiment.foil_separation),...
    num2str(experiment.foil_offset),...
    num2str(experiment.offset_p1),...
    num2str(experiment.offset_h1),...
    num2str(experiment.offset_p2),...
    num2str(experiment.offset_h2),...
    num2str(experiment.Temperature),...
    num2str(experiment.pitch_axis),...
    num2str(experiment.piv_var),...
    num2str(experiment.filt_var),...
    experiment.expf_name,...
    experiment.save_lrs};

% generate prompt window

prompt = {'Enter chord size (in meters): ','Enter span (in meters): ','Enter foil shapes (as string): ', ...
    'Enter Mean wall distance (left, in meters): ','Enter Mean wall distance (right, in meters): ', ...
    'Enter Flume water height(in meters): ', 'Enter anticipated flume frequency (Hz): ', ... 
    'Enter number of foils in experiment: ','Enter foil separation distance (m): ','Enter foil offset distance (m): ', ...
    'Enter Gromit THETA offset from starting position (deg): ','Enter Gromit HEAVE offset from starting position (m): ',...
    'Enter Wallace THETA offset from homing position (deg): ','Enter Wallace HEAVE offset from homing position (m): ',...
    'Enter flume water temperature (from vectrino, in c):   ','Enter foil Pitch Axis','Using PIV? (enter 1 for yes)', ... 
    'Filter data at 60Hz? (0 or 1)','Enter Experiment Folder Name:','Save to LRS? (y/n):'};

name = 'Experiment Configuration';
num_lines = 1; 

% input dialog window
answer = inputdlg(prompt, name, num_lines, defaultanswers);

% store reponses
experiment.chord = str2double(answer{1});
experiment.span = str2double(answer{2});
experiment.foil_shape = char(answer(3));
experiment.Wall_distance_left = str2double(answer{4});
experiment.Wall_distance_right = str2double(answer{5});
experiment.flume_height = str2double(answer{6});
experiment.flume_hertz = str2double(answer{7});
experiment.Number_of_foils = str2double(answer{8});
experiment.foil_separation = str2double(answer{9});
experiment.foil_offset = str2double(answer{10});
experiment.offset_p1 = str2double(answer{11});
experiment.offset_h1 = str2double(answer{12});
experiment.offset_p2 = str2double(answer{13});
experiment.offset_h2 = str2double(answer{14});
experiment.Temperature = str2double(answer{15});
experiment.pitch_axis = str2double(answer{16});
experiment.piv_var = str2double(answer{17});
experiment.filt_var = str2double(answer{18});

% establish save folder
if answer{20} == 'y'
    experiment.fname = ['R:\ENG_Breuer_Shared\group\Joel\main_flume_experimental_data\',answer{19}];
else
    experiment.fname = ['D:\Experiments\',num2str(experiment.Number_of_foils),'foil\',answer{19}];
end

mkdir(experiment.fname);
folder_name = [experiment.fname,'\data'];
mkdir(folder_name); % for the bias measurements

experiment.expf_name = answer{19};

experiment.srate = fs;
experiment.T = 1/fs;

experiment.offset_home = [experiment.offset_p1, experiment.offset_h1, experiment.offset_p2, experiment.offset_h2];

end

