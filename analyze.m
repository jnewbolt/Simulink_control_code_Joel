% This script will run analysis on the data that is the folder/ specified by the variable "filename"

% Load data, see Libraries/Analysis/dataLocations.m for more data storage location strings
datadir = 'C:\Users\Joel\Documents\Globus-BreuerLab@Home\';
trialdir = 'FoilAndCircCyl_D=24,2cm\data\'; namepart1 = '20230524_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
thcknss = 0.054; %cross-stream diameter in meters % Is this necessary? Try to remove

% Type of analysis requested 
singletrial_analysis = 1;
manytrial_analysis = 0;
varyphase = 1;

if manytrial_analysis==1
fstarvector = (0.1:0.02:0.6);%(0.05:0.01:0.14);
Astarvector = (0.0:0.05:1.1);
ftrials = length(fstarvector); Atrials = length(Astarvector);
    if varyphase==1
    fstarvector = 0.09;
    phase12vector = (-180:20:180);
    ftrials = length(phase12vector); 
    end
elseif singletrial_analysis==1
fstarvector = 0.1;
Astarvector = 0;
ftrials = 1; Atrials = 1;
    if varyphase==1
    phase12vector = 0;
    end
end

% % Combine strings to form filename and load last trial to get some necessary variable values from the trial
trialfilename = [datadir,trialdir,namepart1];
trialfiles = dir([datadir,trialdir]);
load([datadir,trialdir,trialfiles(4).name]); % also should remove the need for this

% Define size of each variable to preallocate memory for arrays
variablesPreallocated;

% Loop through trials with different flow speed U
for ftrial = 1:ftrials

% Loop through subtrials with different heave amplitude A
for Atrial = 1:Atrials
    
    loadTrialData;

    % Extract measured quantities
%     [time_star,heave_commanded,heave_measured,heave_star_measured,pitch_measured,force_D,force_L,inertialload_y,...
%     flowspeed_measured,heave_velo,heave_accel] = extract_measurements(transientcycs,freq,T,Prof_out_angle,out);

    % Define timesteps for each subtrial, excluding ramp up/down
    timesteps = length(out);
    timestep_start = round(transient_cycs/(freq*T))+1;
    timestep_end = round(timesteps-transient_cycs/(freq*T));
    times = T*(1:timestep_end-timestep_start+1);
    time_star = times*freq;
    phase12(ftrial,Atrial) = phase;

%     heave_commanded = Prof_out_angle(timestep_start:timestep_end,6);
    heave_commanded = outh2.Data;
    heave_measured = out(timestep_start:timestep_end,4);
    heave_star_measured = heave_measured/thcknss;
    pitch_measured = out(timestep_start:timestep_end,3);
    force_x0 = out(timestep_start:timestep_end,7);
    force_y0 = out(timestep_start:timestep_end,8);
    force_z0 = out(timestep_start:timestep_end,9);
    torque_x0 = out(timestep_start:timestep_end,10);
    torque_y0 = out(timestep_start:timestep_end,11);
    torque_z0 = out(timestep_start:timestep_end,12);
    force_D = force_x0.*cos(pitch_measured) - force_y0.*sin(pitch_measured);
    force_L = force_y0.*cos(pitch_measured) + force_x0.*sin(pitch_measured);
    inertialload_y = out(timestep_start:timestep_end,23);
    flowspeed_measured = abs(out(timestep_start:timestep_end,13));

    heave_velo = movmean((1/T)*gradient(squeeze(heave_measured)),100);
    heave_accel = movmean((1/T)*gradient(squeeze(heave_velo)),100);

    flowspeed_measured_mean(ftrial,Atrial) = mean(flowspeed_measured);
    flowspeed_measured_stdev(ftrial,Atrial) = std(flowspeed_measured);
    Ustar_measured_mean(ftrial,Atrial) = flowspeed_measured_mean(ftrial,Atrial)/(thcknss*freq);

    f_star_commanded(ftrial,Atrial) = thcknss*freq/flowspeed_measured_mean(ftrial,Atrial);
    A_star_commanded(ftrial,Atrial) = (max(heave_commanded)-min(heave_commanded))/(2*thcknss);
    A_star_measured(ftrial,Atrial) = (max(heave_measured)-min(heave_measured))/(2*thcknss);

% Filter force data
    force_L_corrected = force_L+inertialload_y;
    [b,a] = butter(6,10*freq*(2*T),'low'); % butterworth filter 6th order with cut-off frequency at 10*freq
    force_L_corrected_filtered = filtfilt(b,a,squeeze(force_L_corrected)); 
    force_D_filtered = filtfilt(b,a,squeeze(force_D));
    torque_x0_filtered = filtfilt(b,a,squeeze(torque_x0));
    force_scale(ftrial,Atrial) = 0.5*1000*thcknss*foil.span*flowspeed_measured_mean(ftrial,Atrial)^2;
    liftcoef = force_L_corrected_filtered/force_scale(ftrial,Atrial);
    dragcoef = force_D_filtered/force_scale(ftrial,Atrial);
    dragtorquecoef = -torque_x0_filtered/(force_scale(ftrial,Atrial)*foil.span/2);

% Calculate power
    power_scale(ftrial,Atrial) = 0.5*1000*thcknss*foil.span*flowspeed_measured_mean(ftrial,Atrial)^3;
    power_fluid = force_L_corrected_filtered .*heave_velo;
    power_damping = 0*heave_velo.^2;
    power_total = power_fluid + power_damping;
    power_mean(ftrial,Atrial) = mean(power_fluid);
    powercoef = power_fluid/power_scale(ftrial,Atrial);
    powercoef_mean(ftrial,Atrial) = power_mean(ftrial,Atrial)/power_scale(ftrial,Atrial);

% % heave spectrum
    duration = max(time_star/freq)/T;
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(0*window_duration);
    heave_star_hilbert =  hilbert(heave_star_measured);
    [heave_powerspec,f_heave] = pwelch(heave_star_hilbert,window_duration,overlap,[],1/T);
    [max_power,max_index] = max(10*log10(heave_powerspec));
    f_heave_dom = f_heave(max_index);


 % force corrected+filtered spectrum using Welch's method
    duration = round(max(time_star/freq)/T);
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(window_duration*7/8);
    force_hilbert =  hilbert(liftcoef);
    [force_powerspec,f_force] = pwelch(force_hilbert,window_duration,overlap,[],1/T);
    [max_force_power,max_force_index] = max(10*log10(force_powerspec));
    f_force_dom = f_force(max_force_index);  
%     spacing = (f_force(2)-f_force(1))/f;
%     findpeaks(10*log10(force_powerspec),1/spacing,'MinPeakHeight',max_power-20)

% %     % Find phase delay between heave and force
% %     maxlag = round(1/(2*freq*T));  % Maximum lag to calculate xcorr (in timesteps)
% %     [corrs,lags] = xcorr(heave_measured,force_y_corrected_filtered,maxlag);
% % %     delay = freq*T*360*finddelay(heave_measured,force_y_corrected_filtered,maxlag)
% %     hold on
% %     plot(lags*freq*T*360,corrs)
% %     [forcedelay_peak,forcedelay_peaklocs] = findpeaks(corrs,1);
% %     forcedelay_deg = forcedelay_peaklocs*T*freq*360-180; % Delay in degrees
% %     stem(forcedelay_deg,forcedelay_peak);
% %     xlabel('Lag (degrees)'); ylabel('Correlation')
% %     xlim([-180 180])
% 
% %     % Hilbert transform for instantaneous freq
% %     window = floor(duration/10);
% %     noverlap = floor(window/4);
% %     fftpoints = window;
% %     spectrogram(liftcoef,window,noverlap,fftpoints,1/T,'yaxis');
% %     caxis([-30 10])
% %     ylim([0 10*freq])
% %     xlim([0 max(times)])
% 
%    % Find peaks
% spacing = (f_force(2)-f_force(1))/freq;
%     [forcespec_peakpowers,forcespec_peaklocs] = findpeaks(10*log10(force_powerspec),1/spacing, ...
%     'MinPeakHeight',max_force_power-20,'MinPeakDistance',0.5);
%     num_forcespecpeaks = size(forcespec_peaklocs,1);
%     
% % Plot power spectrum
% St = 0.21; % estimated Strouhal number
% f_vortex = St*(flowspeed_measured_mean(ftrial,Atrial)/chord);
% plot_PrescribedMotionPowerSpectrum_MATLABin
% %(freq,f_vortex,f_heave_dom,f_heave,heave_powerspec,f_force,force_powerspec)

% Plot force and motion
    if singletrial_analysis==1
        close all
    figure
    plot_PrescribedMotionForceAndVelocity_MATLABin(time_star,heave_star_measured,heave_velo,liftcoef,...
        dragcoef,power_fluid,num_cyc,dragtorquecoef);
    drawnow
    end
end



end

