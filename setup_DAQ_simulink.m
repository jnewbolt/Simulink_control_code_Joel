%% Setup DAQ - Simulink control
% Run this code before running the experiment in order to:
% -Define user parameters
% -Take tare measurement for force sensors
% -Find zero y-Force angles, aka find zero pitch
% -Take loaded tare measurement (with steady flow on object)

setupEmail() % sets up email to send out when experiment is finished

%% General setup

srate = 1000;
T = 1/srate;

experiment = setupPrompt(srate);
foil = foils_database(experiment.foil_shape);

% expected time delay between Gromit and Wallace (Gromit leading motion)
experiment.motion_delay = 13;
disp(['NOTE: Expected time delay between Gromit and Wallace motions (Gromit leading the motion) is set to ',num2str(experiment.motion_delay),' ms']);

clearvars -except experiment foil srate T

%% Unloaded bias measurement

% check for required parameters:
if ~exist('experiment','var')
    error('Missing necessary variables from workspace')
end

disp(['Finding unloaded bias. Ensure flume is OFF and motors are ON.',newline,'Press any key to continue.'])
pause()

run("find_bias_simulink.m")

clearvars -except experiment foil srate T out raw_encoders raw_wallace raw_gromit bias

%% Find zero pitch

align_ans = input(['Run "find_zero_pitch" for Gromit? y/n + Enter',newline],"s");
if strcmp(align_ans,'y')
    % align gromit
    traverse = 'g';
    disp('Ensure flume is at speed and Gromit is ON. Press any key to continue')
    pause

    % check for required parameters:
    if ~exist('experiment','var') || ~exist('bias','var') || ~exist('traverse','var')
        error('Missing necessary variables from workspace. Vars "experiment", "bias", and "traverse" must be available.')
    end
    
    run("find_zero_pitch_simulink.m")
    
    clearvars -except experiment foil srate T out raw_encoders raw_wallace raw_gromit bias
end

align_ans = input(['Run "find_zero_pitch" for Wallace? y/n + Enter',newline],"s");
if strcmp(align_ans,'y')
    % align wallace
    traverse = 'w';
    disp('Ensure flume is at speed and Wallace is ON Press any key to continue')
    pause
    
    run("find_zero_pitch_simulink.m")
end

% Assign aligned pitch bias to unloaded bias variable
disp('Updating "bias_unloaded" with alignment.')
% bias.pitch(1) = -0.0013; % temporary 20230713
% bias.pitch(2) = 0.7192; % temporary 20230418
bias_unloaded = bias;
clearvars -except experiment foil srate T out raw_encoders raw_wallace raw_gromit bias bias_unloaded

%% Loaded bias measurement

% check for required parameters:
if ~exist('experiment','var') || ~exist('bias_unloaded','var') || ~exist('bias','var')
    % in this case, bias and bias_unloaded are the same, but "find_bias" uses the variable "bias" to run
    error('Missing necessary variables from workspace. Vars "experiment", "bias", and "bias_unloaded" must be available.')
end

disp(['Finding loaded bias. Ensure flume is ON and motors are ON. Bring flowspeed up to experiment value.',newline,'Press any key to continue.'])
pause
run("find_bias_simulink.m")

% Assign bias to loaded bias variable
bias_loaded = bias;
clearvars -except experiment foil srate T out raw_encoders raw_wallace raw_gromit bias bias_unloaded bias_loaded

%% Ready

disp('Done initializing experimental setup.')

% b = timeseries(1);
% set_param('test_model','SimulationCommand','start');

