%% Setup DAQ - Simulink control
% Run this code before running the experiment in order to:
% -Define user parameters
% -Take tare measurement for force sensors
% -Find zero y-Force angles, aka find zero pitch
% -Take loaded tare measurement (with steady flow on object)
%% General setup
setup_email() % sets up email to send out when experiment is finished
ExperimentParameters = setup_prompt(); % prompts user to input parameter values for the experiment

% expected time delay between Gromit and Wallace (Gromit leading motion)
ExperimentParameters.motionDelay = 13;
disp(['NOTE: Expected time delay between Gromit and Wallace motions (Gromit leading the motion) is set to ', ...
    num2str(ExperimentParameters.motionDelay),' ms']);
%% Unloaded bias measurement
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['Finding unloaded bias. Ensure flume is OFF and motors are ON.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
% This is not done with a function call so that the Simulink model can access workspace variables
find_bias_simulink
clearvars -except ExperimentParameters Measurements Biases

%% Find zero pitch
answer = input(['Run "find_zero_pitch" for Gromit? y/n + Enter', newline],"s");
switch answer
    case 'y'
    traverse = 'Gromit';
    disp('Ensure flume is at speed and Gromit is ON. Press any key to continue')
    pause
    % check for required parameters:
    if ~exist('ExperimentParameters','var') || ~exist('Biases','var') || ~exist('Measurements','var')
        error('Missing necessary variables from workspace. Structures "ExperimentParameters", "Biases", and "Measurements" must be available.')
    end
% Run find_zero_pitch_simulink.m script, then clear temporary variables from the workspace
% This is not done with a function call so that the Simulink model can access workspace variables
    find_zero_pitch_simulink
    clearvars -except ExperimentParameters Biases Measurements

end

answer = input(['Run "find_zero_pitch" for Wallace?? y/n + Enter', newline],"s");
switch answer
    case 'y'
    traverse = 'Wallace';
    disp('Ensure flume is at speed and Wallace is ON Press any key to continue')
    pause
    % check for required parameters:
    if ~exist('ExperimentParameters','var') || ~exist('Biases','var') || ~exist('Measurements','var')
        error('Missing necessary variables from workspace. Structures "ExperimentParameters", "Biases", and "Measurements" must be available.')
    end
% Run find_zero_pitch_simulink.m script, then clear temporary variables from the workspace
% This is not done with a function call so that the Simulink model can access workspace variables
    find_zero_pitch_simulink
    clearvars -except ExperimentParameters Biases Measurements
end
    clearvars -except ExperimentParameters Biases Measurements
%% Loaded bias measurement
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['Finding loaded bias. Ensure flume is ON at speed and motors are ON.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
% This is not done with a function call so that the Simulink model can access workspace variables
find_bias_simulink
clearvars -except ExperimentParameters Measurements Biases

%% Ready
disp('Done initializing experimental setup.')
