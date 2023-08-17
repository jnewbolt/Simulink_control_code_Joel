%% Setup DAQ - Simulink control
% Run this code before running the experiment in order to:
% -Define user parameters
% -Take tare measurement for force sensors
% -Find zero y-Force angles, aka find zero pitch
% -Take loaded tare measurement (with steady flow on object)
%% General setup
setup_email() % sets up email to send out when experiment is finished
Parameters = setup_prompt(); % prompts user to input parameter values for the experiment

%% Move motors to starting position
disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure they have clearance then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,Parameters.firstFoilPitchOffsetDegrees); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,Parameters.firstFoilHeaveOffsetMeters); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(0,Parameters.secondFoilPitchOffsetDegrees); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(0,Parameters.secondFoilHeaveOffsetMeters); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Unloaded bias measurement
disp(['Finding unloaded bias. Ensure flume is OFF and motors are ON.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Measurements Biases endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Find zero pitch
answer = input(['Run "find_zero_pitch" for Gromit? y/n + Enter', newline],"s");
switch answer
    case 'y'
    traverse = 'Gromit';
    disp('Ensure flume is at speed and Gromit is ON. Press any key to continue')
    pause
% Run find_zero_pitch_simulink.m script, then clear temporary variables from the workspace
    run('find_zero_pitch_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
    clearvars -except Parameters Biases Measurements  endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

end

answer = input(['Run "find_zero_pitch" for Wallace?? y/n + Enter', newline],"s");
switch answer
    case 'y'
    traverse = 'Wallace';
    disp('Ensure flume is at speed and Wallace is ON Press any key to continue')
    pause
    % Run find_zero_pitch_simulink.m script, then clear temporary variables from the workspace
    run('find_zero_pitch_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
    clearvars -except Parameters Biases Measurements endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW
end
    clearvars -except Parameters Biases Measurements endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW
    
%% Move to new start positions based on find zero pitch
disp(['The traverses will be moved to their starting positions.',newline, ...
    'Make sure they have clearance then press any key to continue'])
pause()
% Run the move to center of the flume
[startPitchDegG, endPitchDegG] = deal(0,Parameters.firstFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0); %#ok<ASGLU> 
[startPitchDegW, endPitchDegW] = deal(endPitchDegW,Parameters.secondFoilZeroPitchDegrees); %#ok<ASGLU> 
[startHeaveMetersW, endHeaveMetersW] = deal(endHeaveMetersW,endHeaveMetersW); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW
%% Redo the bias measurement at the new start positions
disp(['Repeating unloaded bias. Ensure flume is OFF and motors are ON.',newline,'Press any key to continue.'])
pause()
% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Biases endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Loaded bias measurement
disp(['Finding loaded bias. Ensure flume is ON at speed and motors are ON.',newline,'Press any key to continue.'])
pause()
% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Biases endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Move motors to home position
disp('The traverses will move to their home positions.')
% Run the move to start positions
[startPitchDegG, endPitchDegG] = deal(endPitchDegG,0); 
[startHeaveMetersG, endHeaveMetersG] = deal(0,0);
[startPitchDegW, endPitchDegW] = deal(endPitchDegW,0); 
[startHeaveMetersW, endHeaveMetersW] = deal(endPitchDegW,0);
run('move_to_position')
clearvars -except Parameters Biases endPitchDegG endHeaveMetersG endPitchDegW endHeaveMetersW

%% Ready
disp('Done initializing experimental setup.')
