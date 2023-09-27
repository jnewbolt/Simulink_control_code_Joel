function foils = foils_database(foilType1,foilType2, print_properties)
% This function acts as a database for the different foil pairs that may be
% used in the Flume. It requires a string of the desired foil/foils as an
% input.
%
% v2.0.0, Handy, 03/17/2023
% -------------------------------------------------------------------------
% Inputs:
%
% foiltype = string containing the two or three digit code for the desired
%            foil.
%
% print_properties = string, "y". If defined this function will print the
%                    properties of the selected foil to the command line.
%
% -------------------------------------------------------------------------
% Outputs:
%
% foil = structure containing properties of the selected foil:
%           > chord (in m), AR (Aspect Ratio), mass1 (in kg), mass2 (in
%             case there is a foil pair), span (in m), profile (shape),
%             material.
% -------------------------------------------------------------------------
% Examples:
%
% foil = foils_database('A1')
%
% foil = foils_database('A1', [])
%
% foil = foils_database('A1', 'y')
%
% -------------------------------------------------------------------------
% NOTES:
% -------------------------------------------------------------------------

% Properties of each foil: chord (m), thickness (m), span (m), mass (kg), shape, material
FoilProps.A2 = struct('chord',0.0762,'thickness',NaN,'span',0.457,'mass',0.928,'shape','Rectangular','material','Aluminum'); % Yunxing's foils
FoilProps.AE3 = struct('chord',0.061,'thickness',NaN,'span',0.366,'mass',0.538,'shape','Rectangular','material','Aluminum');% Eric's main medium foils
FoilProps.C1 = struct('chord',0.054,'thickness',0.054,'span',0.447,'mass',0.386,'shape','Cylindrical','material','Carbon fiber');% Joel's cylinder
FoilProps.EC1 = struct('chord',0.0594,'thickness',0.0238,'span',0.401,'mass',0.302,'shape','Elliptical cylinder','material','PLA');% Joel's elliptical cylinder
FoilProps.V1 = struct('chord',0.0535,'thickness',0.0265,'span',0.401,'mass',0.306,'shape','VibrissaBeem50x','material','PLA');% Joel's vibrissae model
FoilProps.C2 = struct('chord',0.0538,'thickness',0.0538,'span',0.400,'mass',0.598,'shape','Cylindrical','material','ABS');% Joel's cylinder Xometry
FoilProps.None = 'No 2nd foil';

foils.firstFoil = FoilProps.(foilType1);
foils.secondFoil = FoilProps.(foilType2);

if exist('print_properties','var')
    if strcmp(print_properties,'y')
        foil % prints out the data of the selected foil
    elseif ~isempty(print_properties)
        error("Wrong argument for input var 'print_properties', must either be 'y' or empty.")
    end
end
% Old table of foil properties
% FoilProperties =...
%       [" Chord [m]: ";" AspectRatio :";" Mass1 [kg]  :";" Mass2 [kg]  :";" Profile     :";   " Material    :"  ]; % COMMENTS:
% E1 =  [     0.1;              3.5;           0.450;           0.450;         "Elliptical";      "Carbon fiber"    ]; % Built by someone
% A1 =  [     0.1;              3.5;           1.026;           1.174;         "Rectangular";     "Aluminum"        ]; % Nick's foils
% A2 =  [   3*0.0254;           6.0;           0.928;           0.924;         "Rectangular";     "Aluminum"        ]; % Yunxing's foils
% C1 =  [     0.054;            8.28;          0.386;           0.0;           "Cylindrical";     "Carbon fiber"    ]; % Joel's cylinder
% V1 =  [     0.0535;           7.50;          0.306;           0.0;           "VibrissaeBeem50x";"PLA3DprintWepoxy"]; % Joel's vibrissae model
% EC1 = [     0.0594;           6.81;          0.302;           0.0;           "Elliptical";      "PLA3DprintWepoxy"]; % Joel's vibrissae model
% A3E = [     0.061;            6.0;           0.538;           0.538;         "Rectangular";     "Aluminum"        ]; % Eric's main medium foils
% F1E = [   3*0.0254;           6.0;           0.614;           0.0;           "NACA0022";        "Silicone"        ]; % Ian's fishtail foil with shaft extension
%
% foils = table(FoilProperties, E1, A1, A2, C1, V1, EC1, A3E, F1E); % constructing a table out of the foil data
% 
% selected_data = foils.(foiltype); % identifies the column of the selected foils
% 
% foil.chord = str2double(selected_data(1)); % asigns first column's value to the chord variable
% foil.AR    = str2double(selected_data(2)); % asigns second column's value to the aspect ratio var
% foil.mass1 = str2double(selected_data(3)); % asigns third column's value to the leading foil's mass var
% foil.mass2 = str2double(selected_data(4)); % asigns fourth column's value to the trailing foil's mass var
% foil.span  = foil.chord*foil.AR; % calculates the span
% foil.profile  = selected_data(5); % asigns fifth column's value to the geometric profile of the foils
% foil.material = selected_data(6); % asigns sixth column's value to the geometric profile of the foils
end
