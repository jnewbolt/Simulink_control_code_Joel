function    [time_star,heave_commanded,heave_measured,heave_star_measured,pitch_measured,force_D,force_L,inertialload_y,...
    flowspeed_measured,heave_velo,heave_accel] = extract_measurements(transientcycs,freq,T,Prof_out_angle,out)
    % Define timesteps for each subtrial, excluding ramp up/down
    chord = 0.0594;
    timesteps = length(out);
    timestep_start = round(transientcycs/(freq*T))+1;
    timestep_end = round(timesteps-transientcycs/(freq*T));
    times = T*(1:timestep_end-timestep_start+1);
    time_star = times*freq;

    heave_commanded = Prof_out_angle(timestep_start:timestep_end,6);
    heave_measured = out(timestep_start:timestep_end,6);
    heave_star_measured = heave_measured/chord;
    pitch_measured = out(timestep_start:timestep_end,5);
    force_x0 = out(timestep_start:timestep_end,7);
    force_y0 = out(timestep_start:timestep_end,8);
    force_D = force_y0.*cos(pitch_measured) - force_x0.*sin(pitch_measured);
    force_L = force_x0.*cos(pitch_measured) + force_y0.*sin(pitch_measured);
    inertialload_y = out(timestep_start:timestep_end,23);
    flowspeed_measured = abs(out(timestep_start:timestep_end,13));

    heave_velo = movmean((1/T)*gradient(squeeze(heave_measured)),100);
    heave_accel = movmean((1/T)*gradient(squeeze(heave_velo)),100);

end