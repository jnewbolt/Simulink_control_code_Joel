function [] = plot_power_spectrum(flowspeedMetersPerSecMean,chordMeters,freq,freqSpec,powerSpec,freqCutoff,plotTitle)

   figure('Position', [100 80 1600 900])
    hold on; grid on;



    % Plot heave spectrum

%     plot([f_heave_dom f_heave_dom],[-60 30],'Color','red','DisplayName','Frequency dominant');
%     plot(f_heave,10*log10(heave_powerspec));
%     plot([f_force_dom f_force_dom],[-60 30],'Color','red','LineWidth',4);
    plot([1 1],[-60 60],'Color','red','LineWidth',4,'DisplayName','Flapping frequency');
    plot([freqCutoff/freq freqCutoff/freq],[-60 60],'Color','blue','LineWidth',4,'DisplayName','Filter cutoff');

%     plot([flume_hertz/freq,flume_hertz/freq],[-60,60],'color','cyan','LineStyle','--','LineWidth',4)
    plot(freqSpec/freq,10*log10(powerSpec),'Color','black','LineWidth',4,'DisplayName','PSD'); % plot force spectrum
%     scatter(forcespec_peaklocs,forcespec_peakpowers,'Marker','o');

    % Freq. of predicted vortex shedding
%     St=0.2;
%     fVortex = St*flowspeedMetersPerSecMean/chordMeters;
%     plot([fVortex/freq fVortex/freq],[-60 60],'color','blue','LineStyle','--','LineWidth',4,'DisplayName','Predicted vortex freq.')

    xlabel('Frequency (f/f_{prescribed})')
    ylabel('Force PSD (dB/Hz)')
    xlim([0 20]) 
    ylim([-80 20])
    xticks(linspace(0,20,21))
    title(plotTitle)
    legend()
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off

end