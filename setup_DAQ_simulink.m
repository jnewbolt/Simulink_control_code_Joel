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
    'Make sure MOTORS ARE ON with CLEARANCE to move, then press any key to continue.'])
pause()
% Run the move to center of the flume
[startPitchDegG, MotorPositions.endPitchDegG] = deal(0,Parameters.pitchOffsetDegG); %#ok<ASGLU> 
[startHeaveMetersG, MotorPositions.endHeaveMetersG] = deal(0,Parameters.heaveOffsetMetersG); %#ok<ASGLU> 
[startPitchDegW, MotorPositions.endPitchDegW] = deal(0,Parameters.pitchOffsetDegW); %#ok<ASGLU> 
[startHeaveMetersW, MotorPositions.endHeaveMetersW] = deal(0,Parameters.heaveOffsetMetersW); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters MotorPositions

%% Unloaded bias measurement
disp(['Finding UNLOADED BIAS. Ensure flume is OFF.',newline,'Press any key to continue.'])
pause()

% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Measurements Biases MotorPositions

%% Find zero pitch
answer = input(['Run "find_zero_pitch" for Gromit? y/n + Enter', newline],"s");
switch answer
    case 'y'
    traverse = 'Gromit';
    disp('Turn FLUME ON. Then press any key to continue')
    pause
% Run find_zero_pitch_simulink.m script, then clear temporary variables from the workspace
    run('find_zero_pitch_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
    clearvars -except Parameters Biases Measurements MotorPositions

end

answer = input(['Run "find_zero_pitch" for Wallace?? y/n + Enter', newline],"s");
switch answer
    case 'y'
    traverse = 'Wallace';
    disp('Ensure FLUME is ON. Then press any key to continue')
    pause
    % Run find_zero_pitch_simulink.m script, then clear temporary variables from the workspace
    run('find_zero_pitch_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
    clearvars -except Parameters Biases Measurements MotorPositions
end
    clearvars -except Parameters Biases Measurements MotorPositions
    
%% Move to new start positions based on find zero pitch
disp('The traverses will be moved to their zero-pitch positions.')
% Run the move to center of the flume
[startPitchDegG, MotorPositions.endPitchDegG] = deal(MotorPositions.endPitchDegG,Parameters.pitchOffsetDegG); %#ok<ASGLU> 
[startHeaveMetersG, MotorPositions.endHeaveMetersG] = deal(MotorPositions.endHeaveMetersG,MotorPositions.endHeaveMetersG); %#ok<ASGLU> 
[startPitchDegW, MotorPositions.endPitchDegW] = deal(MotorPositions.endPitchDegW,Parameters.pitchOffsetDegW); %#ok<ASGLU> 
[startHeaveMetersW, MotorPositions.endHeaveMetersW] = deal(MotorPositions.endHeaveMetersW,MotorPositions.endHeaveMetersW); %#ok<ASGLU> 
run('move_to_position') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters MotorPositions
%% Redo the bias measurement at the new start positions
disp(['Repeating UNLOADED BIAS. Turn FLUME OFF.',newline,'Press any key to continue.'])
pause()
% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Biases MotorPositions

%% Loaded bias measurement
disp(['Finding LOADED BIAS. Turn FLUME ON.',newline,'Press any key to continue.'])
pause()
% Run find_bias_simulink.m script, then clear temporary variables from the workspace
run('find_bias_simulink') % This is not done with a function call so that the Simulink model can access workspace variables
clearvars -except Parameters Biases MotorPositions

%% Move motors to home position
disp('The traverses will move to their home positions.')
% Run the move to start positions
[startPitchDegG, MotorPositions.endPitchDegG] = deal(MotorPositions.endPitchDegG,0); 
[startHeaveMetersG, MotorPositions.endHeaveMetersG] = deal(MotorPositions.endHeaveMetersG,0);
[startPitchDegW, MotorPositions.endPitchDegW] = deal(MotorPositions.endPitchDegW,0); 
[startHeaveMetersW, MotorPositions.endHeaveMetersW] = deal(MotorPositions.endHeaveMetersW,0);
run('move_to_position')
clearvars -except Parameters Biases MotorPositions

%% Ready
disp('Done with experiment setup!')
