% This script will run analysis on the data that is the folder/ specified by the variable "filename"

%% Load data, see Libraries/Analysis/dataLocations.m for more data storage location strings
datadir = 'D:\Experiments\Data\FoilAndVib4';%
%% Some alternate important data locations:
% datadir = 'C:\Users\Joel\Documents\Globus-BreuerLab@Home\FoilAndVib_D=75cm';
% datadir = 'R:\ENG_Breuer_Shared\group\JoelNewbolt\ExperimentalData\FlexFoilAndCyl_25-Sep-2023_16-47-28';
% trialdir = 'FoilAndVib_D=24,2cm\data\'; namepart1 = '20230525_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
% thcknss = 0.0265;  %cross-stream diameter in meters % Is this necessary? Try to remove
% trialdir = 'FoilAndCircCyl_D=24,2cm\data\'; namepart1 = '20230524_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
% thcknss = 0.054; %cross-stream diameter in meters % Is this necessary? Try to remove
% trialdir = 'FoilSymmetryTest2\data\'; namepart1 = '20230616_PrescribedMotion_p1='; namepart2 = 'deg_h1='; namepart3 = 'cm_ph=0deg';
% thcknss = 0.054; chord_foil = 0.075; span_foil =0.45;%cross-stream diameter in meters % Is this necessary? Try to remove
% trialdir = 'FoilSymmetryTestSmallFoil\data\'; namepart1 = '20230620_PrescribedMotion_p1='; namepart2 = 'deg_h1='; namepart3 = 'cm_ph=0deg';
% thcknss = 0.054; %chord_foil = 0.06; span_foil =0.365;%cross-stream diameter in meters % Is this necessary? Try to remove
% trialdir = 'vib50xBeem\data\'; namepart1 = 'vib_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'CircCyl_20220919\data\'; namepart1 = 'CylPowerMap_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'EllipticalCyl_04-Jul-2022_16_7_4\data\'; namepart1 = 'EllipticalCyl_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'VibManyFreq_27-Oct-2022_18_52_34\data\'; namepart1 = 'Vib_pitch=0deg,f='; namepart2='Hz,A=';
% filename = 'Data\20220620_foilandvib\vary_phase12\foilandvib_pitch=0deg,f='; name2='Hz,A='; name3='cm,phase12=';
% datadir = 'C:\Users\Joel\Documents\Globus-BreuerLab@Home\';
% trialdir = 'FoilAndVib_close\data\'; namepart1 = '20230427_PrescribedMotion_p2=0deg_h2='; namepart2 = 'c_ph='; namepart3 = 'deg';
%%
%% Type of analysis requested 
singleTrialAnalysis = 1;
manyTrialAnalysis = 0;
varyphase = 0;
varypitch1 = 0;

createplotForcesVsTimeG = 0;
createplotForcesVsTimeW = 1;
createplotForcesVsDisplacementsG = 0;
createplotForcesVsDisplacementsW = 0;
createplotForceSpectrumG = 0;
createplotForceSpectrumW = 1;
createplotEnergyMap = 0;
createGIF = 0;
createVideo = 0;

if manyTrialAnalysis == 0
    firstTrial = 1;
    nTrials = firstTrial;
end

%% Find the data files
dir_fullname = [datadir,'\data\'];
trialfiles = dir([dir_fullname,'*.mat']);
trialfiles=trialfiles(contains({trialfiles.name},'trial')); 
trialfiles = natsortfiles(trialfiles);
if manyTrialAnalysis == 1
firstTrial = 1;
nTrials = length(trialfiles);
end
% 
% % % % Combine strings to form filename and load last trial to get some necessary variable values from the trial
% trialfilename = [datadir,trialdir,namepart1];
% trialfiles = dir([datadir,trialdir]);
% load([datadir,trialdir,trialfiles(4).name]); % also should remove the need for this

%%


% if manytrial_analysis==1
% fstarvector = (0:2:10);%(0.05:0.01:0.14);
% Astarvector = (0:0.1:0.5);%*(0.0265/0.0535);
% ftrials = length(fstarvector); Atrials = length(Astarvector);
%     if varyphase==1
%     fstarvector = (0:2:10);
%     phase12vector = (0:2:10);
%     ftrials = length(phase12vector); 
%     end
% elseif singletrial_analysis==1
% fstarvector = 10;
% Astarvector = 0.5*(0.0265/0.0535);
% ftrials = 1; Atrials = 1;
%     if varyphase==1
%     phase12vector = (-180:20:180); 
%     ftrials = length(phase12vector); 
%     end
% end


%% Define the necessary variables in order to pre-allocate memory
freqStarCmdW = nan(nTrials,1);
heaveAmpStarCmdW = nan(nTrials,1);
heaveAmpStarMeasuredW = nan(nTrials,1);
phase12 = nan(nTrials,1);
flowspeedMetersPerSecMean = nan(nTrials,1);
flowspeedMetersPerSecStddev = nan(nTrials,1);
UstarMean = nan(nTrials,1);
powerMeanG = nan(nTrials,1);
powerCoefMeanG = nan(nTrials,1);
powerScale = nan(nTrials,1);
flowspeed_measured_p2p = nan(nTrials,1);
f_force_filtered_dom = nan(nTrials,1);
forceScaleW = nan(nTrials,1);
forceScaleG = nan(nTrials,1);
delay = nan(nTrials,1);
num_forcespecpeaks = nan(nTrials,1);
liftcoef_alltrials_G = nan(nTrials,1);
dragcoef_alltrials_G = nan(nTrials,1);
inertialload_alltrials = nan(nTrials,1);
torqueLDscaleG = nan(nTrials,1);
torquezScaleG = nan(nTrials,1);
torqueLDscaleW = nan(nTrials,1);
torquezScaleW = nan(nTrials,1);
angleOfAttackMaxDegW = nan(nTrials,1);
videoIndex = 1;
%% Loop through trials

for iTrial = firstTrial:nTrials

%% Load the requested trial data

    try
        trialname=[dir_fullname,trialfiles(iTrial).name];
        load(trialname);
    catch
        disp(['Failed to load ',trialname])
    end

%% Extract measured quantities
    P = Parameters; M = Measurements;
    % Define timesteps for each subtrial, excluding ramp up/down
    nTimesteps = length(Measurements.pitchDegreesG);
    timestepFirst = round(nTransientCycsG*P.sampleRate/freqG)+1;
    timestepLast = round(nTimesteps-nTransientCycsG*P.sampleRate/freqG);
    times = (1:timestepLast-timestepFirst+1)/P.sampleRate;
    timeStarG = times*freqG;
    timeStarW = times*freqW;
    phase12(iTrial) = phaseLagWbehindG;

    chordMetersG = Parameters.Foils.foilG.chord;
    chordMetersW = Parameters.Foils.foilW.chord;
    spanMetersG = Parameters.Foils.foilG.span;
    spanMetersW = Parameters.Foils.foilW.span;
    heaveMetersCropG = M.heaveMetersG(timestepFirst:timestepLast); 
    heaveMetersCropW = M.heaveMetersW(timestepFirst:timestepLast); 
    heaveStarG = heaveMetersCropG/chordMetersG;
    heaveStarW = heaveMetersCropW/chordMetersW;
    pitchRadsCropG = (pi/180)*M.pitchDegreesG(timestepFirst:timestepLast); 
    pitchRadsCropW = (pi/180)*M.pitchDegreesW(timestepFirst:timestepLast); 
    forcexW = M.forcesNewtonsW(timestepFirst:timestepLast,1);
    forceyW = M.forcesNewtonsW(timestepFirst:timestepLast,2);
    forcezW = M.forcesNewtonsW(timestepFirst:timestepLast,3);
    torquexW = M.forcesNewtonsW(timestepFirst:timestepLast,4);
    torqueyW = M.forcesNewtonsW(timestepFirst:timestepLast,5);
    torquezW = M.forcesNewtonsW(timestepFirst:timestepLast,6);
    forcexG = M.forcesNewtonsG(timestepFirst:timestepLast,1);
    forceyG = M.forcesNewtonsG(timestepFirst:timestepLast,2);
    forcezG = M.forcesNewtonsG(timestepFirst:timestepLast,3);
    torquexG = M.forcesNewtonsG(timestepFirst:timestepLast,4);
    torqueyG = M.forcesNewtonsG(timestepFirst:timestepLast,5);
    torquezG = M.forcesNewtonsG(timestepFirst:timestepLast,6);
    inertiaLoadyFiltG = M.forceInertialLoadG(timestepFirst:timestepLast,1);
    flowSpeedMetersPerSec = M.flowMetersPerSecondVectrino(timestepFirst:timestepLast,1);

    % Flow speed statistical quantities
    flowspeedMetersPerSecMean(iTrial) = abs(mean(flowSpeedMetersPerSec));
    flowspeedMetersPerSecStddev(iTrial) = std(flowSpeedMetersPerSec);
    UstarMean(iTrial) = flowspeedMetersPerSecMean(iTrial)/(chordMetersG*freqG);

    % Forces in flow coordinate system
    forceDragG = forcexG.*cos(pitchRadsCropG) - forceyG.*sin(pitchRadsCropG);
    forceLiftG = forceyG.*cos(pitchRadsCropG) + forcexG.*sin(pitchRadsCropG);    
    forceDragW = forcexW.*cos(pitchRadsCropW) - forceyW.*sin(pitchRadsCropW);
    forceLiftW = forceyW.*cos(pitchRadsCropW) + forcexW.*sin(pitchRadsCropW);
    torqueDragG = torqueyG.*cos(pitchRadsCropG) + torquexG.*sin(pitchRadsCropG);
    torqueLiftG = -torquexG.*cos(pitchRadsCropG) + torqueyG.*sin(pitchRadsCropG);
    torqueDragW = torqueyW.*cos(pitchRadsCropW) + torquexW.*sin(pitchRadsCropW);
    torqueLiftW = -torquexW.*cos(pitchRadsCropW) + torqueyW.*sin(pitchRadsCropW);
    % Velocities and accelerations from position derivatives
    movemeanPoints = 20;
    pitchVelocityDegPerSecG = movmean(gradient(squeeze(pitchRadsCropG)*P.sampleRate),movemeanPoints);
    pitchAccelDegPerSecSqG = movmean(gradient(squeeze(pitchVelocityDegPerSecG)*P.sampleRate),movemeanPoints);
    heaveVelocityMetersPerSecG = movmean(gradient(squeeze(heaveMetersCropG)*P.sampleRate),movemeanPoints);
    heaveAccelMetersPerSecSqG = movmean(gradient(squeeze(heaveVelocityMetersPerSecG)*P.sampleRate),movemeanPoints);
    pitchVelocityDegPerSecW = movmean(gradient(squeeze(pitchRadsCropW)*P.sampleRate),movemeanPoints);
    pitchAccelDegPerSecSqW = movmean(gradient(squeeze(pitchVelocityDegPerSecW)*P.sampleRate),movemeanPoints);
    heaveVelocityMetersPerSecW = movmean(gradient(squeeze(heaveMetersCropW)*P.sampleRate),movemeanPoints);
    heaveAccelMetersPerSecSqW = movmean(gradient(squeeze(heaveVelocityMetersPerSecW)*P.sampleRate),movemeanPoints);
    % Dimensionless quantities
    freqStarCmdW(iTrial) = chordMetersW*freqG/flowspeedMetersPerSecMean(iTrial);
    heaveAmpStarCmdW(iTrial) = (max(trajHeaveMetersW)-min(trajHeaveMetersW))/(2*chordMetersW);
    heaveAmpStarMeasuredW(iTrial) = (max(heaveMetersCropW)-min(heaveMetersCropW))/(2*chordMetersW);
    angleOfAttackMaxDegW(iTrial) = (180/pi)*atan2(2*pi*freqW*heaveAmpMetersW,flowspeedMetersPerSecMean(iTrial)) - pitchAmpDegW;
%% Filter force data
    forceLiftGcorrected = forceLiftG;%+inertialload_y_G;
    freqCutoff = 10*freqG;
    [b,a] = butter(6,freqCutoff*2/P.sampleRate,'low'); % butterworth filter 6th order with cut-off frequency at 10*freq
    forceLiftFiltG = filtfilt(b,a,squeeze(forceLiftGcorrected));
    forceLiftFiltW = filtfilt(b,a,squeeze(forceLiftW)); 
    forceDragFiltG = filtfilt(b,a,squeeze(forceDragG));
    forceDragFiltW = filtfilt(b,a,squeeze(forceDragW));
    torqueDragFiltG = filtfilt(b,a,squeeze(torqueDragG));
    torqueDragFiltW = filtfilt(b,a,squeeze(torqueDragW));
    torqueLiftFiltG = filtfilt(b,a,squeeze(torqueLiftG));
    torqueLiftFiltW = filtfilt(b,a,squeeze(torqueLiftW));
    torquezFiltG = filtfilt(b,a,squeeze(torquezG));
    torquezFiltW = filtfilt(b,a,squeeze(torquezW));
    inertiaLoadyFiltG = filtfilt(b,a,squeeze(inertiaLoadyFiltG));
    %% Nondimnesionalize forces and torques using characteristic scales
    forceScaleG(iTrial) = 0.5*1000*chordMetersG*spanMetersG*flowspeedMetersPerSecMean(iTrial)^2; 
    forceScaleW(iTrial) = 0.5*1000*chordMetersW*spanMetersW*flowspeedMetersPerSecMean(iTrial)^2;
    torqueLDscaleG(iTrial) = forceScaleG(iTrial)*spanMetersG;
    torquezScaleG(iTrial) = forceScaleG(iTrial)*chordMetersG;
    torqueLDscaleW(iTrial) = forceScaleW(iTrial)*spanMetersW;
    torquezScaleW(iTrial) = forceScaleW(iTrial)*chordMetersW;
    liftCoefG = forceLiftFiltG/forceScaleG(iTrial);
    liftCoefW = forceLiftFiltW/forceScaleW(iTrial);
    dragCoefG = forceDragFiltG/forceScaleG(iTrial);
    dragCoefW = forceDragFiltW/forceScaleW(iTrial);
    torqueLiftCoefG = torqueLiftFiltG/torqueLDscaleG(iTrial);
    torqueDragCoefG = torqueDragFiltG/torqueLDscaleG(iTrial);
    torquezCoefG = torquezFiltG/torquezScaleG(iTrial);
    torqueLiftCoefW = torqueLiftFiltW/torqueLDscaleW(iTrial);
    torqueDragCoefW = torqueDragFiltW/torqueLDscaleW(iTrial);
    torquezCoefW = torquezFiltW/torquezScaleW(iTrial);

%% Calculate power
    powerScale(iTrial) = 0.5*1000*chordMetersG*spanMetersG*flowspeedMetersPerSecMean(iTrial)^3;
    powerFluid = forceLiftFiltG .*heaveVelocityMetersPerSecG;
    powerDamping = 0*heaveVelocityMetersPerSecG.^2;
    powerTotal = powerFluid + powerDamping;
    powerMeanG(iTrial) = mean(powerFluid);
    powerCoefG = powerFluid/powerScale(iTrial);
    powerCoefMeanG(iTrial) = powerMeanG(iTrial)/powerScale(iTrial);
    
%% Spectrum analysis
    %Gromit Heave spectrum
    durationG = max(timeStarG/freqG)*P.sampleRate;
    windowDurationG = round(durationG/2); % size of hilbert windows measured in samples
    overlapG = round(0*windowDurationG);
    heaveHilbertG =  hilbert(heaveStarG);
    [heaveSpecPowerG,heaveSpecFreqG] = pwelch(heaveHilbertG,windowDurationG,overlapG,[],P.sampleRate);
    [maxHeavePowerG,maxHeavePowerIndexG] = max(10*log10(heaveSpecPowerG));
    freqHeaveDominantG = heaveSpecFreqG(maxHeavePowerIndexG);

    % Wallace Heave spectrum
    durationW = max(timeStarW/freqW)*P.sampleRate;
    windowDurationW = round(durationW/2); % size of hilbert windows measured in samples
    overlapW = round(0*windowDurationW);
    heaveStarHilbertW =  hilbert(heaveStarW);
    [heaveSpecPowerW,heaveSpecFreqW] = pwelch(heaveStarHilbertW,windowDurationW,overlapW,[],P.sampleRate);
    [maxHeavePowerW,maxHeavePowerIndexW] = max(10*log10(heaveSpecPowerW));
    freqHeaveDominantW = heaveSpecFreqW(maxHeavePowerIndexW);

    %Gromit Pitch spectrum
    durationG = max(timeStarG/freqG)*P.sampleRate;
    windowDurationG = round(durationG/2); % size of hilbert windows measured in samples
    overlapG = round(0*windowDurationG);
    pitchHilbertG =  hilbert(pitchRadsCropG);
    [pitchSpecPowerG,pitchSpecFreqG] = pwelch(pitchHilbertG,windowDurationG,overlapG,[],P.sampleRate);
    [maxPitchPowerG,maxPitchPowerIndexG] = max(10*log10(pitchSpecPowerG));
    freqPitchDominantG = pitchSpecFreqG(maxPitchPowerIndexG);

    %Wallace pitch spectrum
    durationW = max(timeStarW/freqW)*P.sampleRate;
    windowDurationW = round(durationW/2); % size of hilbert windows measured in samples
    overlapW = round(0*windowDurationW);
    pitchHilbertW =  hilbert(pitchRadsCropW);
    [pitchSpecPowerW,pitchSpecFreqW] = pwelch(pitchHilbertW,windowDurationW,overlapW,[],P.sampleRate);
    [maxPitchPowerW,maxPitchPowerIndexW] = max(10*log10(pitchSpecPowerW));
    freqPitchDominantW = pitchSpecFreqW(maxPitchPowerIndexW);

 % force corrected+filtered spectrum using Welch's method
    durationG = round(max(timeStarG/freqG)*P.sampleRate);
    windowDurationG = round(durationG/2); % size of hilbert windows measured in samples
    overlapG = round(windowDurationG*1/2);
    forceHilbertG =  hilbert(liftCoefG);%hilbert(liftcoef_G);
    [forcePowerSpecG,freqForceSpecG] = pwelch(forceHilbertG,hanning(windowDurationG),overlapG,[],P.sampleRate);
    [maxForcePowerG,maxForceIndexG] = max(10*log10(forcePowerSpecG));
    freqForceDominantG = freqForceSpecG(maxForceIndexG);  
 % force corrected+filtered spectrum using Welch's method
    durationW = round(max(timeStarW/freqW)*P.sampleRate);
    windowDurationW = round(durationW/2); % size of hilbert windows measured in samples
    overlapW = round(windowDurationW*1/2);
    forceHilbertW =  hilbert(liftCoefW);%hilbert(liftcoef_W);
    [forcePowerSpecW,freqForceSpecW] = pwelch(forceHilbertW,hanning(windowDurationW),overlapW,[],P.sampleRate);
    [maxForcePowerW,maxForceIndexW] = max(10*log10(forcePowerSpecW));
    freqForceDominantW = freqForceSpecW(maxForceIndexW);  
    
    %% Phase average quantities
    samplesPerCycW = round(P.sampleRate/freqW);
    timeStarPhaseAvgW = (1/samplesPerCycW:1/samplesPerCycW:1);
    pitchRadsPhaseAvgW = nan(samplesPerCycW,1);
    pitchRadsPhaseStddevW = nan(samplesPerCycW,1);
    heaveStarPhaseAvgW = nan(samplesPerCycW,1);
    heaveStarPhaseStddevW = nan(samplesPerCycW,1);
    liftCoefPhaseAvgW = nan(samplesPerCycW,1);
    liftCoefPhaseStddevW = nan(samplesPerCycW,1);
    dragCoefPhaseAvgW =nan(samplesPerCycW,1);
    dragCoefStddevW = nan(samplesPerCycW,1);
        for iSample = 1:samplesPerCycW
        pitchRadsPhaseAvgW(iSample) = mean(pitchRadsCropW(iSample:samplesPerCycW:end));
        pitchRadsPhaseStddevW(iSample) = std(pitchRadsCropW(iSample:samplesPerCycW:end));
        heaveStarPhaseAvgW(iSample) = mean(heaveStarW(iSample:samplesPerCycW:end));
        heaveStarPhaseStddevW(iSample) = std(heaveStarW(iSample:samplesPerCycW:end));
        liftCoefPhaseAvgW(iSample) = mean(liftCoefW(iSample:samplesPerCycW:end));
        liftCoefPhaseStddevW(iSample) = std(liftCoefW(iSample:samplesPerCycW:end));
        dragCoefPhaseAvgW(iSample) = mean(dragCoefW(iSample:samplesPerCycW:end));
        dragCoefStddevW(iSample) = std(dragCoefW(iSample:samplesPerCycW:end));
        end
    samplesPerCycG = round(P.sampleRate/freqG);
    timeStarPhaseAvgG = (0:freqG/P.sampleRate:1);
    pitchRadsPhaseAvgG = nan(samplesPerCycG,1);
    pitchRadsPhaseStddevG = nan(samplesPerCycG,1);
    heaveStarPhaseAvgG = nan(samplesPerCycG,1);
    heaveStarPhaseStddevG = nan(samplesPerCycG,1);
    liftCoefPhaseAvgG = nan(samplesPerCycG,1);
    liftCoefPhaseStddevG = nan(samplesPerCycG,1);
    dragCoefPhaseAvgG =nan(samplesPerCycG,1);
    dragCoefStddevG = nan(samplesPerCycG,1);
        for iSample = 1:samplesPerCycG
        pitchRadsPhaseAvgG(iSample) = mean(pitchRadsCropG(iSample:samplesPerCycG:end));
        pitchRadsPhaseStddevG(iSample) = std(pitchRadsCropG(iSample:samplesPerCycG:end));
        heaveStarPhaseAvgG(iSample) = mean(heaveStarG(iSample:samplesPerCycG:end));
        heaveStarPhaseStddevG(iSample) = std(heaveStarG(iSample:samplesPerCycG:end));
        liftCoefPhaseAvgG(iSample) = mean(liftCoefG(iSample:samplesPerCycG:end));
        liftCoefPhaseStddevG(iSample) = std(liftCoefG(iSample:samplesPerCycG:end));
        dragCoefPhaseAvgG(iSample) = mean(dragCoefG(iSample:samplesPerCycG:end));
        dragCoefStddevG(iSample) = std(dragCoefG(iSample:samplesPerCycG:end));
        end
    %%
%     spacing = (f_force(2)-f_force(1))/f;
%     findpeaks(10*log10(force_powerspec1/spacing,'MinPeakHeight',max_power-20)

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
if createplotForceSpectrumG == 1
    plotTitle = ['Object behind foil phase diff ',num2str(phaseLagWbehindG,2),'deg and +/-',num2str(heaveAmpMetersG/chordMetersG,2),'chord'];
    % plot_PrescribedMotionPowerSpectrum_MATLABin
    plot_power_spectrum(flowspeedMetersPerSecMean(iTrial),chordMetersG,freqG,freqForceSpecG,forcePowerSpecG,plotTitle)
end
if createplotForceSpectrumW == 1
    plotTitle = ['Foil w/ heave amplitude ',num2str(heaveAmpMetersW/chordMetersW,2),'{\it c}',newline,...
        'and pitch amplitude ',num2str(pitchAmpDegW,2),'deg',newline,...
        'AoA max = ',num2str(angleOfAttackMaxDegW(iTrial),2),' deg'];% plot_PrescribedMotionPowerSpectrum_MATLABin
    plot_power_spectrum(flowspeedMetersPerSecMean(iTrial),chordMetersW,freqW,freqForceSpecW,forcePowerSpecW,plotTitle)
end
%% Plot force and motion
% if (mod(iTrial,1)==0 && iTrial>6*length(phaseLagWbehindGvec)+1 && ...
%          iTrial<7*length(phaseLagWbehindGvec)+1 && singleTrialAnalysis==1) || ...
%          createGIF == 0 % only plot every Nth figure

    if createplotForcesVsTimeG == 1
% Plot lift and drag coefficients and heave position for Gromit
    titlePlots = ['Object downstream',newline,'w/ heave amplitude ',num2str(heaveAmpMetersG/chordMetersG,2),'{\it c} and phase diff ',num2str(phaseLagWbehindG,2),'deg'];
%     plotForceTorqueDisplacementVsTime(timeStar,pitchRadsCropG,heaveStarG,liftCoefG,...
%         dragCoefG,powerFluid,torqueLiftCoefG,torqueDragCoefG,torquezCoefG,nCyclesG,...
%         titlePlots);
    plot_forces_vs_time_bounded(timeStarPhaseAvgG,pitchRadsPhaseAvgG,pitchRadsPhaseStddevG,...
        heaveStarPhaseAvgG,heaveStarPhaseStddevG,...
        liftCoefPhaseAvgG,liftCoefPhaseStddevG,dragCoefPhaseAvgG,dragCoefStddevG,...
        powerFluid,nCyclesG,titlePlots);
    end
    if createplotForcesVsTimeW==1
% Plot lift and drag coefficients and heave position for Wallace
    titlePlots = ['Foil w/ heave amplitude ',num2str(heaveAmpMetersW/chordMetersW,2),'{\it c}',newline,...
        'and pitch amplitude ',num2str(pitchAmpDegW,2),'deg',newline,...
        'AoA max = ',num2str(angleOfAttackMaxDegW(iTrial),2),' deg'];
%     plotForceTorqueDisplacementVsTime(timeStar,pitchRadsCropW,heaveStarW,liftCoefW,...
%         dragCoefW,powerFluid,torqueLiftCoefW,torqueDragCoefW,torquezCoefW,nCyclesG,...
%         titlePlots);
    plot_forces_vs_time_bounded(timeStarPhaseAvgW,pitchRadsPhaseAvgW,pitchRadsPhaseStddevW,...
        heaveStarPhaseAvgW,heaveStarPhaseStddevW,...
        liftCoefPhaseAvgW,liftCoefPhaseStddevW,dragCoefPhaseAvgW,dragCoefStddevW,...
        powerFluid,nCyclesG,titlePlots);
    end
% Plot lift and drag coefficients vs heave and angular position for Gromit
    if createplotForcesVsDisplacementsG==1
        titlePlots = ['Object behind foil +/-',num2str(pitchAmpDegG,2),'deg and +/-',num2str(heaveAmpMetersG/chordMetersG,2),'{\it c}'];
        plotForceTorqueVsDisplacement(timeStarG,P.sampleRate,freqG,pitchRadsCropG,heaveStarG,liftCoefG,...
            torquezCoefG,titlePlots)
    end
    if createplotForcesVsDisplacementsW==1
        titlePlots = ['Foil ahead of object +/-',num2str(pitchAmpDegW,2),'deg and +/-',num2str(heaveAmpMetersW/chordMetersW,2),'{\it c}'];
        plotForceTorqueVsDisplacement(timeStarG,P.sampleRate,freqG,pitchRadsCropW,heaveStarW,liftCoefW,...
            torquezCoefW,titlePlots)
    end

    if createVideo == 1
    % Make video frames from set of figures that are currently being generated
    drawnow
    fig = gcf;
    frame = getframe(fig);
    im{videoIndex} = frame2im(frame);
    videoIndex = videoIndex+1;
    end
    if manyTrialAnalysis == 1 % close figures after loading for many trial analysis
%     close all
    end

% end
end
%% Save .mp4 video of plots
if createVideo == 1
    filenameVideo = "testVideo.mp4";
    nImages = length(im);
    v = VideoWriter(filenameVideo,'MPEG-4');
    v.Quality = 95;
    v.FrameRate = 2;
    frameHoldTimeSec = 2;
    open(v);
    for idx = 1:nImages
            for repeatFrame = 1:frameHoldTimeSec*v.FrameRate
            frame = im{idx};
            writeVideo(v,frame)
            end
    end
    close(v);
end
%% Plot Power coefficient contour
% plot_PowerCoefContour_vsphasediffandampl

% Plot vs phase
% f_star_sorted = f_star_commanded;
% U_star_sorted = 1./f_star_sorted;
% Repeat the first row because data is periodic -180deg = 180deg

% phase12_sorted = [phase12;-phase12(1,:)];
% A_star_sorted = [A_star_measured;A_star_measured(1,:)];
% powercoef_mean_sorted = [powercoef_mean;powercoef_mean(1,:)];


if createplotEnergyMap ==1
% Plot acceleration limit
acc_limit = 4.9; % Acceleration limit in m/s^2
v_limit = 0.5;

figure
hold on
% Plot Cp vs. A* and phase12
powerCoefMeanReshape = reshape(powerCoefMeanG,length(phaseLagWbehindGvec),length(heaveAmpMetersGvec));
contourf(phaseLagWbehindGvec,heaveAmpMetersGvec,powerCoefMeanReshape',120,'LineStyle','none')
clim([-1 0.2])
colormap(bluewhitered)

contour(phaseLagWbehindGvec,heaveAmpMetersGvec,powerCoefMeanReshape',[-1e-6 -1e-6],'LineWidth',4,'LineColor','k','LineStyle','--')
% scatter(phase_vec,H2star_vec,60,'.','k')
grid on
xlabel('Phase between foil and object (degrees)')
ylabel('{\it A} * = {\it A/d}')
xlim([-190 180])
ylim([-0.04 1.12])
set(gca, 'Layer', 'top')
grid off

c=colorbar();
c.Label.String = '{\it C}_P';
% c.Label.Interpreter = 'Latex';
plotTitle = ['Object behind foil, foil pitch +/-',num2str(pitchAmpDegG,2),'deg and +/-',num2str(heaveAmpMetersG/chordMetersG,2),'chord'];
title(plotTitle);
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off
end