function out = convert_output(raw_encoders, raw_force_wallace, raw_force_gromit, raw_vectrino, raw_accelmeter, ref_signal, bias, range_x, offset_home)

% % Required inputs:
% raw_encoders
% raw_wallace
% raw_gromit

% % Optional inputs:
% bias
% range_x
% offset_p2
% offset_h2

%% Encoders

if ~exist('range_x', 'var') || isempty(range_x)
    range_x = 1:size(raw_encoders,1);
end

if ~exist('bias','var') || isempty(bias)
    bias.Wallace = 0;
    bias.Gromit = 0;
end

if ~exist('offset_home','var') || isempty(offset_home)
    offset_p2 = 0;
    offset_h2 = 0;
else
    offset_p2 = offset_home(3);
    offset_h2 = offset_home(4);
end

raw_encoders(raw_encoders(:,1:2)>1e6) = raw_encoders(raw_encoders(:,1:2)>1e6)-2^32;

out(:,1) = (raw_encoders(range_x,1)).*2*pi/10000;
out(:,2) = -raw_encoders(range_x,2).*0.0254/8000;
out(:,3) = deg2rad(conv_encodertheta_traverse(raw_encoders(range_x,3), 0)) - deg2rad(offset_p2);
out(:,4) = conv_encodery_traverse(raw_encoders(range_x,4), 0) - offset_h2;

out(:,5) = ref_signal(range_x);

%% Force transducers

out(:,7:12) = conv_output_wallace_force(raw_force_wallace(range_x,:),bias.Wallace);
out(:,17:22) = conv_output_gromit_force(raw_force_gromit(range_x,:),bias.Gromit);

%% Vectrino

out(:,13:16) = (raw_vectrino(range_x,:)-2.5)*2/5;

%% Accelerometer

out(:,23) = (0.6+0.306)*9.81*(raw_accelmeter(range_x,:)-bias.accmeter); % Mass of mount + mass of vibrissa*(m/s/s/Volt)*accmeter_voltage

end

