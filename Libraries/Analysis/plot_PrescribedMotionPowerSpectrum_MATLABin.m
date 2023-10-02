
%         % Plot heave spectrum
%     figure
%     hold on; grid on;
%     plot([f_heave_dom f_heave_dom],[-60 30],'Color','red');
%     plot(f_heave,10*log10(heave_powerspec));
% %     plot([f_force_dom f_force_dom],[-60 30],'Color','red','DisplayName','Frequency dominant');
% %     plot(f_force,10*log10(force_powerspec),'Color','black','DisplayName','Force spectrum'); % plot force spectrum
%     xlabel('Frequency (Hz)')
%     ylabel('Heave PSD (dB/Hz)')
%     xlim([0 10]) 
%     ylim([0 30])
% %     title(['CPS ring down, simulated: m = ',num2str(M),' kg, f_{natural} = ',num2str(f_nat),' Hz, \zeta = ',num2str(zeta)]);
%     hold off

   figure('Position', [100 80 1600 900])
    hold on; grid on;

    % Plot heave spectrum

%     plot([f_heave_dom f_heave_dom],[-60 30],'Color','red','DisplayName','Frequency dominant');
%     plot(f_heave,10*log10(heave_powerspec));
%     plot([f_force_dom f_force_dom],[-60 30],'Color','red','LineWidth',4);
    plot([1 1],[-60 60],'Color','blue','LineStyle','--','LineWidth',4);

%     % freq of flapping foil
%     plot([1.02/freq 1.02/freq],[-60 60],'Color','red','LineStyle','--');
%     % foil harmonics
%     plot([2*1.02/freq 2*1.02/freq],[-60 60],'Color','red','LineStyle',':');
%     plot([3*1.02/freq 3*1.02/freq],[-60 60],'Color','red','LineStyle',':');

    % Freq. of predicted vortex shedding
    St=0.2;
    fVortex = St*flowspeedMetersPerSecMean/chordMetersG;
    plot([fVortex/freqG fVortex/freqG],[-60 60],'color','red','LineStyle','--','LineWidth',4)

%     plot([flume_hertz/freq,flume_hertz/freq],[-60,60],'color','cyan','LineStyle','--','LineWidth',4)
    plot(freqForceSpec/freqG,10*log10(forcePowerSpec),'Color','black','LineWidth',4); % plot force spectrum
%     scatter(forcespec_peaklocs,forcespec_peakpowers,'Marker','o');
    xlabel('Frequency (f/f_{prescribed})')
    ylabel('Force PSD (dB/Hz)')
    xlim([0 20]) 
    ylim([-80 20])
    xticks(linspace(0,10,11))
    title(plotTitle)
%     title(['CPS ring down, simulated: m = ',num2str(M),' kg, f_{natural} = ',num2str(f_nat),' Hz, \zeta = ',num2str(zeta)]);
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off

