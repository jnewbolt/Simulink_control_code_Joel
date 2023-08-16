%% Setup DAQ - Simulink control
% Run this code before running the experiment in order to:
% -Define user parameters
% -Take tare measurement for force sensors
% -Find zero y-Force angles, aka find zero pitch
% -Take loaded tare measurement (with steady flow on object)
%% General setup
setup_email() % sets up email to send out when experiment is finished
ExperimentParameters = setup_prompt(); % prompts user to input parameter values for the experiment

%% Move motors to starting position
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure they have clearance then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,ExperimentParameters.firstFoilPitchOffsetDegrees); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,ExperimentParameters.firstFoilHeaveOffsetMeters); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(0,ExperimentParameters.secondFoilPitchOffsetDegrees); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(0,ExperimentParameters.secondFoilHeaveOffsetMeters); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except ExperimentParameters
%% Unloaded bias measurement
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['Finding unloaded bias. Ensure flume is OFF and motors are ON.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
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
    run('find_zero_pitch_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
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
    run('find_zero_pitch_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
    clearvars -except ExperimentParameters Biases Measurements
end
    clearvars -except ExperimentParameters Biases Measurements
%% Move to new start positions based on find zero pitch
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure they have clearance then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,ExperimentParameters.firstFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(endPitchDegW,ExperimentParameters.secondFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(endHeaveMetersW,endHeaveMetersW); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except ExperimentParameters
%% Redo the bias measurement at the new start positions
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['Repeating unloaded bias. Ensure flume is OFF and motors are ON.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except ExperimentParameters Measurements Biases

%% Loaded bias measurement
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp(['Finding loaded bias. Ensure flume is ON at speed and motors are ON.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except ExperimentParameters Measurements Biases

%% Move motors to home position
% check for required parameters:
if ~exist('ExperimentParameters','var')
    error('Missing necessary variables from workspace')
end

disp('The traverses will move to their home positions.')
% Run the move to start positions
[startPitchDegG, endPitchDegG] = deal(endPitchDegG,0); 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0);
[startPitchDegW, endPitchDegW] = deal(endPitchDegW,0); 
[startHeaveMetersW, endHeaveMetersW] = deal(endPitchDegW,0);
run('move_to_position')
clearvars -except ExperimentParameters Measurements Biases

%% Ready
disp('Done initializing experimental setup.')
