function [t_ramp, ramp_p, ramp_h] = ramp_fn(ramp_time, T, bias, offset_home, traverse)
    % used by both traverses in order to ramp up and down from starting position to offset and bias
    switch traverse
        case 'g'
            offset_p = offset_home(1);
            offset_h = offset_home(2);
            bias_pitch = bias.pitch(1);
            bias_heave = bias.heave(1);
        case 'w'
            offset_p = offset_home(3);
            offset_h = offset_home(4);
            bias_pitch = bias.pitch(2);
            bias_heave = bias.heave(2);
    end

    Npts_ramp = ramp_time/T;
    t_ramp = T*(0:Npts_ramp-1); % dimensional ramping time vector
    
    ramp_p = (offset_p+bias_pitch)*(0.5*(1-cos( pi*(0:Npts_ramp-1)/Npts_ramp)))';
    ramp_h = (offset_h+bias_heave)*(0.5*(1-cos( pi*(0:Npts_ramp-1)/Npts_ramp)))';

end