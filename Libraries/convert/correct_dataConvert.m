
Gromit_Cal = [-0.021960016340,0.004190267995,-0.186118155718,-17.483510971069,1.224225640297,16.644166946411;
    0.965337395668,19.793544769287,-0.233169600368,-10.133704185486,-0.765193104744,-9.669131278992;
    29.811275482178,0.075847007334,30.349569320679,0.450678527355,29.886123657227,0.429086029530;
    0.018749127164,0.455327808857,-1.043430089951,-0.245171427727,1.023705482483,-0.209242567420;
    1.182301402092,-0.001399876550,-0.586799740791,0.394854426384,-0.641057908535,-0.38928997516;
    -0.028681563213,-0.638982295990,-0.011541210115,-0.659547865391,-0.044132679701,-0.627585828304];
Gromit_Cal_inv = inv(Gromit_Cal);
Wallace_Cal = [-0.358090000000000,  -0.0161800000000000, 0.338480000000000,-17.4630700000000,    0.890470000000000,  16.8772400000000;
           0.155930000000000,  21.5597800000000,   -0.785970000000000,-10.0453800000000,    0.430010000000000,  -9.68643000000000;
          29.8995100000000,     0.375370000000000, 29.6925600000000,    0.457840000000000, 30.5462600000000,     0.601450000000000;
           0.00323000000000000, 0.492950000000000, -1.03800000000000,  -0.243410000000000,  1.06263000000000,   -0.203400000000000;
           1.21724000000000,    0.0112800000000000,-0.589940000000000,  0.391800000000000, -0.631210000000000,  -0.395370000000000;
          -0.00657000000000000,-0.716070000000000,  0.0176200000000000,-0.663460000000000, -0.0425100000000000, -0.644210000000000];
Wallace_Cal_inv = inv(Wallace_Cal);

% Correct data conversion for bias trials
trialfiles = dir('R:\ENG_Breuer_Shared\jnewbolt\DAQandMotorControl\Data\FoilAndVib_close\data\bias*.mat');

for i=2:size(trialfiles,1)
    
    load([trialfiles(i).folder,'\',trialfiles(i).name]);
    out(:,1) = out(:,1)*(360/(2*pi)); % convert from incorrect units to radians

    out(:,7:12) = out(:,7:12)*Gromit_Cal_inv'*Wallace_Cal';
    out(:,17:22) = out(:,17:22)*Wallace_Cal_inv'*Gromit_Cal';
    save([trialfiles(i).folder,'\',trialfiles(i).name]);
end

% % Correct data conversion for experiment trials
trialfiles = dir('R:\ENG_Breuer_Shared\jnewbolt\DAQandMotorControl\Data\FoilAndVib_close\data\2023*.mat');

for i=1:size(trialfiles,1)
    
    load([trialfiles(i).folder,'\',trialfiles(i).name]);
    out(:,1) = out(:,1)*(360/(2*pi)); % convert from incorrect units to radians

    out(:,7:12) = out(:,7:12)*Gromit_Cal_inv'*Wallace_Cal';
    out(:,17:22) = out(:,17:22)*Wallace_Cal_inv'*Gromit_Cal';
    save([trialfiles(i).folder,'\',trialfiles(i).name]);
end
