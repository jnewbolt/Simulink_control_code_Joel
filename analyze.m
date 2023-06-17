% This script will run analysis on the data that is the folder/ specified by the variable "filename"

%% Load data, see Libraries/Analysis/dataLocations.m for more data storage location strings
datadir = 'C:\Users\Joel\Documents\Globus-BreuerLab@Home\';
trialdir = 'FoilAndVib_D=24,2cm\data\'; namepart1 = '20230525_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
thcknss = 0.0265;  %cross-stream diameter in meters % Is this necessary? Try to remove
trialdir = 'FoilAndCircCyl_D=24,2cm\data\'; namepart1 = '20230524_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
thcknss = 0.054; %cross-stream diameter in meters % Is this necessary? Try to remove
trialdir = 'FoilSymmetryTest2\data\'; namepart1 = '20230616_PrescribedMotion_p1='; namepart2 = 'deg_h1='; namepart3 = 'cm_ph=0deg';
thcknss = 0.054; chord_foil = 0.075; span_foil =0.45;%cross-stream diameter in meters % Is this necessary? Try to remove
%% Some alternate important data locations:
% trialdir = 'vib50xBeem\data\'; namepart1 = 'vib_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'CircCyl_20220919\data\'; namepart1 = 'CylPowerMap_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'EllipticalCyl_04-Jul-2022_16_7_4\data\'; namepart1 = 'EllipticalCyl_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'VibManyFreq_27-Oct-2022_18_52_34\data\'; namepart1 = 'Vib_pitch=0deg,f='; namepart2='Hz,A=';
% filename = 'Data\20220620_foilandvib\vary_phase12\foilandvib_pitch=0deg,f='; name2='Hz,A='; name3='cm,phase12=';
% datadir = 'R:\ENG_Breuer_Shared\jnewbolt\DAQandMotorControl\Data\';
% trialdir = 'FoilAndVib_close\data\'; namepart1 = '20230427_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
%%
%% Type of analysis requested 
singletrial_analysis = 1;
manytrial_analysis = 1;
varyphase = 1;
varypitch1 = 1;

if manytrial_analysis==1
fstarvector = (0:2:10);%(0.05:0.01:0.14);
Astarvector = (0:0.1:0.5);%*(0.0265/0.0535);
ftrials = length(fstarvector); Atrials = length(Astarvector);
    if varyphase==1
    fstarvector = 0.09;
    phase12vector = (0:2:10);
    ftrials = length(phase12vector); 
    end
elseif singletrial_analysis==1
fstarvector = 0.1;
Astarvector = 1.1;%*(0.0265/0.0535);
ftrials = 1; Atrials = 1;
    if varyphase==1
    phase12vector = 0;%(-180:20:180); 
    ftrials = length(phase12vector); 
    end
end

% % Combine strings to form filename and load last trial to get some necessary variable values from the trial
trialfilename = [datadir,trialdir,namepart1];
trialfiles = dir([datadir,trialdir]);
load([datadir,trialdir,trialfiles(4).name]); % also should remove the need for this

%% Define the necessary variables in order to pre-allocate memory
f_star_commanded_W = nan(ftrials,Atrials);
A_star_commanded_W = nan(ftrials,Atrials);
A_star_measured_W = nan(ftrials,Atrials);
phase12 = nan(ftrials,Atrials);
flowspeed_measured_mean = nan(ftrials,Atrials);
flowspeed_measured_stdev = nan(ftrials,Atrials);
power_mean = nan(ftrials,Atrials);
powercoef_mean = nan(ftrials,Atrials);
power_scale = nan(ftrials,Atrials);
flowspeed_measured_p2p = nan(ftrials,Atrials);
f_force_filtered_dom = nan(ftrials,Atrials);
force_scale_W = nan(ftrials,Atrials);
force_scale_G = nan(ftrials,Atrials);
delay = nan(ftrials,Atrials);
num_forcespecpeaks = nan(ftrials,Atrials);
liftcoef_alltrials = nan(ftrials,Atrials);
dragcoef_alltrials = nan(ftrials,Atrials);
inertialload_alltrials = nan(ftrials,Atrials);
torqueLD_scale_G = nan(ftrials,Atrials);
torquez_scale_G = nan(ftrials,Atrials);
torqueLD_scale_W = nan(ftrials,Atrials);
torquez_scale_W = nan(ftrials,Atrials);
%% Loop through trials with different flow speed U
% Loop through subtrials with different heave amplitude A
for Atrial = 1:Atrials
for ftrial = 1:ftrials

%% Load the requested trial data
%     if varyphase==0
%         trialname = [trialfilename,num2str(fvector(ftrial),3),namepart2,num2str(Avector(Atrial),3),'cm.mat'];
%     elseif varyphase==1
%         trialname = [trialfilename,num2str(Astarvector(Atrial),3),namepart2,num2str(phase12vector(ftrial)),'deg.mat'];
%     end
% 
%     try
%         load(trialname,'transient_cycs','out','profs','freq','phase')
%     catch
%         disp(['Failed to load ',trialname])
%     end
    if varyphase==0
        trialname = [trialfilename,num2str(fstarvector(ftrial),3),namepart2,num2str(Astarector(Atrial),3),'cm.mat'];
    elseif varyphase==1
        trialname = [trialfilename,num2str(Astarvector(Atrial),3),namepart2,num2str(phase12vector(ftrial)),'deg.mat'];
    end
    if varypitch1 ==1
        trialname = [trialfilename,num2str(P1star_vec(ftrial)),namepart2,num2str(H1star_vec(Atrial)*chord_foil),namepart3,'.mat'];
    end

    try
        load(trialname,'transient_cycs','out','profs','freq','phase')
    catch
        disp(['Failed to load ',trialname])
    end
%% Extract measured quantities
    % Define timesteps for each subtrial, excluding ramp up/down
    timesteps = length(out);
    timestep_start = round(transient_cycs/(freq*T))+1;
    timestep_end = round(timesteps-transient_cycs/(freq*T));
    times = T*(1:timestep_end-timestep_start+1);
    time_star = times*freq;
    phase12(ftrial,Atrial) = phase;

    heave_commanded_G = outh1.Data;
    heave_commanded_W = outh2.Data;
    pitch_commanded_G = outp1.Data;
    pitch_commanded_W = outp2.Data;
    heave_measured_G = out(timestep_start:timestep_end,2);
    heave_measured_W = out(timestep_start:timestep_end,4);
    heave_star_measured_G = heave_measured_G/chord_foil; % Compared to thin Aluminum foil chord of 7.5cm
    heave_star_measured_W = heave_measured_W/thcknss;
    pitch_measured_G = -1/pi*out(timestep_start:timestep_end,1);
    pitch_measured_W = -1/pi*out(timestep_start:timestep_end,3);
    force_x0_W = out(timestep_start:timestep_end,7);
    force_y0_W = out(timestep_start:timestep_end,8);
    force_z0_W = out(timestep_start:timestep_end,9);
    torque_x0_W = out(timestep_start:timestep_end,10);
    torque_y0_W = out(timestep_start:timestep_end,11);
    torque_z0_W = out(timestep_start:timestep_end,12);
    force_x0_G = out(timestep_start:timestep_end,17);
    force_y0_G = out(timestep_start:timestep_end,18);
    force_z0_G = out(timestep_start:timestep_end,19);
    torque_x0_G = out(timestep_start:timestep_end,20);
    torque_y0_G = out(timestep_start:timestep_end,21);
    torque_z0_G = out(timestep_start:timestep_end,22);
    force_D_G = force_x0_G.*cos(pitch_measured_G) - force_y0_G.*sin(pitch_measured_G);
    force_L_G = force_y0_G.*cos(pitch_measured_G) + force_x0_G.*sin(pitch_measured_G);    
    force_D_W = force_x0_W.*cos(pitch_measured_W) - force_y0_W.*sin(pitch_measured_W);
    force_L_W = force_y0_W.*cos(pitch_measured_W) + force_x0_W.*sin(pitch_measured_W);
    torque_D_G = torque_y0_G.*cos(pitch_measured_G) + torque_x0_G.*sin(pitch_measured_G);
    torque_L_G = -torque_x0_G.*cos(pitch_measured_G) + torque_y0_G.*sin(pitch_measured_G);
    torque_D_W = torque_y0_W.*cos(pitch_measured_W) + torque_x0_W.*sin(pitch_measured_W);
    torque_L_W = -torque_x0_W.*cos(pitch_measured_W) + torque_y0_W.*sin(pitch_measured_W);
    
    inertialload_y_G = out(timestep_start:timestep_end,23);
    flowspeed_measured = abs(out(timestep_start:timestep_end,13));

    pitch_velo_G = movmean((1/T)*gradient(squeeze(pitch_measured_G)),100);
    pitch_accel_G = movmean((1/T)*gradient(squeeze(pitch_velo_G)),100);
    heave_velo_G = movmean((1/T)*gradient(squeeze(heave_measured_G)),100);
    heave_accel_G = movmean((1/T)*gradient(squeeze(heave_velo_G)),100);
    heave_velo_W = movmean((1/T)*gradient(squeeze(heave_measured_W)),100);
    heave_accel_W = movmean((1/T)*gradient(squeeze(heave_velo_W)),100);

    flowspeed_measured_mean(ftrial,Atrial) = mean(flowspeed_measured);
    flowspeed_measured_stdev(ftrial,Atrial) = std(flowspeed_measured);
    Ustar_measured_mean(ftrial,Atrial) = flowspeed_measured_mean(ftrial,Atrial)/(thcknss*freq);

    f_star_commanded_W(ftrial,Atrial) = thcknss*freq/flowspeed_measured_mean(ftrial,Atrial);
    A_star_commanded_W(ftrial,Atrial) = (max(heave_commanded_W)-min(heave_commanded_W))/(2*thcknss);
    A_star_measured_W(ftrial,Atrial) = (max(heave_measured_W)-min(heave_measured_W))/(2*thcknss);
%% Filter force data
    force_L_W_corrected = force_L_W+inertialload_y_G;
    [b,a] = butter(6,10*freq*(2*T),'low'); % butterworth filter 6th order with cut-off frequency at 10*freq
    force_L_G_corrected_filtered = filtfilt(b,a,squeeze(force_L_G));
    force_L_W_corrected_filtered = filtfilt(b,a,squeeze(force_L_W_corrected)); 
    force_D_G_filtered = filtfilt(b,a,squeeze(force_D_G));
    force_D_W_filtered = filtfilt(b,a,squeeze(force_D_W));
    torque_D_G_filtered = filtfilt(b,a,squeeze(torque_D_G));
    torque_D_W_filtered = filtfilt(b,a,squeeze(torque_D_W));
    torque_L_G_filtered = filtfilt(b,a,squeeze(torque_L_G));
    torque_L_W_filtered = filtfilt(b,a,squeeze(torque_L_W));
    torque_z_G_filtered = filtfilt(b,a,squeeze(torque_z0_G));
    torque_z_W_filtered = filtfilt(b,a,squeeze(torque_z0_W));
    %% Nondimnesionalize forces and torques using characteristic scales
    force_scale_G(ftrial,Atrial) = 0.5*1000*chord_foil*span_foil*flowspeed_measured_mean(ftrial,Atrial)^2; % for 7.5cm chord foil with 45 cm span
    force_scale_W(ftrial,Atrial) = 0.5*1000*thcknss*foil.span*flowspeed_measured_mean(ftrial,Atrial)^2;
    torqueLD_scale_G(ftrial,Atrial) = force_scale_G(ftrial,Atrial)*span_foil;
    torquez_scale_G(ftrial,Atrial) = force_scale_G(ftrial,Atrial)*chord_foil;
    torqueLD_scale_W(ftrial,Atrial) = force_scale_W(ftrial,Atrial)*foil.span;
    torquez_scale_W(ftrial,Atrial) = force_scale_W(ftrial,Atrial)*thcknss;
    liftcoef_G = force_L_G_corrected_filtered/force_scale_G(ftrial,Atrial);
    liftcoef_W = force_L_W_corrected_filtered/force_scale_W(ftrial,Atrial);
    dragcoef_G = force_D_G_filtered/force_scale_G(ftrial,Atrial);
    dragcoef_W = force_D_W_filtered/force_scale_W(ftrial,Atrial);
    torqueliftcoef_G = torque_L_G_filtered/torqueLD_scale_G(ftrial,Atrial);
    torquedragcoef_G = torque_D_G_filtered/torqueLD_scale_G(ftrial,Atrial);
    torquezcoef_G = torque_z_G_filtered/torquez_scale_G(ftrial,Atrial);
    torqueliftcoef_W = torque_L_W_filtered/torqueLD_scale_W(ftrial,Atrial);
    torquedragcoef_W = torque_D_W_filtered/torqueLD_scale_W(ftrial,Atrial);
    torquezcoef_W = torque_z_W_filtered/torquez_scale_W(ftrial,Atrial);

%% Calculate power
    power_scale(ftrial,Atrial) = 0.5*1000*thcknss*foil.span*flowspeed_measured_mean(ftrial,Atrial)^3;
    power_fluid = force_L_W_corrected_filtered .*heave_velo_W;
    power_damping = 0*heave_velo_W.^2;
    power_total = power_fluid + power_damping;
    power_mean(ftrial,Atrial) = mean(power_fluid);
    powercoef = power_fluid/power_scale(ftrial,Atrial);
    powercoef_mean(ftrial,Atrial) = power_mean(ftrial,Atrial)/power_scale(ftrial,Atrial);

%% heave spectrum
    duration = max(time_star/freq)/T;
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(0*window_duration);
    heave_star_hilbert =  hilbert(heave_star_measured_W);
    [heave_powerspec,f_heave] = pwelch(heave_star_hilbert,window_duration,overlap,[],1/T);
    [max_power,max_index] = max(10*log10(heave_powerspec));
    f_heave_dom = f_heave(max_index);


 % force corrected+filtered spectrum using Welch's method
    duration = round(max(time_star/freq)/T);
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(window_duration*7/8);
    force_hilbert =  hilbert(liftcoef_W);
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

%% Plot power spectrum
% St = 0.21; % estimated Strouhal number
% f_vortex = St*(flowspeed_measured_mean(ftrial,Atrial)/chord);
% plot_PrescribedMotionPowerSpectrum_MATLABin(freq,f_vortex,f_heave_dom,f_heave,heave_powerspec,f_force,force_powerspec)

%% Plot force and motion
    if singletrial_analysis==1
        % Plot lift and drag coefficients and heave position for Wallace
%     plotForceAndPosition(time_star,heave_star_measured_W,heave_velo,liftcoef_W,...
%         dragcoef_W,power_fluid,num_cyc,dragtorquecoef_W,'Vibrissa in wake'); 
%     disp(['Phase diff ',num2str(phase),' Avg Power coef ', num2str(powercoef_mean(ftrial,Atrial))]);pause(5)
    % Plot lift and drag coefficients and heave position for Gromit
    titlePlots = ['Flapping foil +/-',num2str(P1star_vec(ftrial)),'deg and +/-',num2str(H1star_vec(Atrial)),'chord'];
%     plotForceTorqueDisplacementVsTime(time_star,pitch_measured_G,heave_star_measured_G,liftcoef_G,...
%         dragcoef_G,power_fluid,torqueliftcoef_G,torquedragcoef_G,torquezcoef_G,num_cyc,...
%         titlePlots);

    plotForceTorqueVsDisplacement(time_star,T,freq,pitch_measured_G,heave_star_measured_G,liftcoef_G,...
        torquezcoef_G,titlePlots)
    % Make a gif
    drawnow
    fig = gcf;
    frame = getframe(fig);
%     idx = Atrial+(ftrial-1)*Atrials;
    idx = ftrial+(Atrial-1)*ftrials;
    im{idx} = frame2im(frame);

    
%     plotTorqueAndPosition(time_star,pitch_measured_G,heave_star_measured_G,torqueliftcoef_G,...
%         torquedragcoef_G,torquezcoef_G,power_fluid,num_cyc,titlePlots)
    end
end

end
%% Save .gif of plots
    filenameGIF = "testAnimated.gif"; % Specify the output file name
    nImages = length(im);
    for idx = 1:nImages
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filenameGIF,"gif","LoopCount",Inf,"DelayTime",1);
        else
            imwrite(A,map,filenameGIF,"gif","WriteMode","append","DelayTime",1);
        end
    end
    close all
%% Plot Power coefficient contour
% Plot vs phase
% f_star_sorted = f_star_commanded;
% U_star_sorted = 1./f_star_sorted;
% Repeat the first row because data is periodic -180deg = 180deg
% phase12_sorted = [phase12;-phase12(1,:)];
% A_star_sorted = [A_star_measured;A_star_measured(1,:)];
% powercoef_mean_sorted = [powercoef_mean;powercoef_mean(1,:)];

% Plot acceleration limit
% acc_limit = 4.9; % Acceleration limit in m/s^2
% v_limit = 0.5;
% 
% hold on
% % Plot Cp vs. A* and phase12
% contourf(phase12,A_star_measured_W,powercoef_mean,120,'LineStyle','none')
% caxis([-1 0.6])
% colormap(bluewhitered)
% 
% contour(phase12,A_star_measured_W,powercoef_mean,[-1e-6 -1e-6],'LineWidth',4,'LineColor','k','LineStyle','--')
% scatter(phase12,A_star_measured_W,60,'.','k')
% grid on
% xlabel('Phase between foil and vibrissa (degrees)')
% ylabel('{\it A} * = {\it A/d}')
% xlim([-190 180])
% ylim([-0.04 1.12])
% set(gca, 'Layer', 'top')
% grid off
% 
% c=colorbar();
% c.Label.String = '{\it C}_P';
% % c.Label.Interpreter = 'Latex';
% set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
% 
% hold off