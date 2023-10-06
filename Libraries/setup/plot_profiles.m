%% plot profiles

function plot_profiles(timeVector,pitchDegW,heaveMetersW, pitchDegG, heaveMetersG, refSig)

    figure
    subplot(2,1,1)
    plot(timeVector,pitchDegW,"DisplayName","Pitch 1"); hold on;
    plot(timeVector,pitchDegG,"DisplayName","Pitch 2");
    if exist('data_ref', 'var')
        plot(timeVector,refSig*270,'r--',"DisplayName","Reference");
    end
    hold off;
    ylim([-90,270])
    xlabel('time (seconds)')
    ylabel('degrees')
    legend()

    subplot(2,1,2)
    plot(timeVector,heaveMetersW,"DisplayName","Heave 1"); hold on;
    plot(timeVector,heaveMetersG,"DisplayName","Heave 2");
    if exist('data_ref', 'var')
        plot(timeVector,refSig,'r--',"DisplayName","Reference");
    end
    hold off;
    ylim([-0.15,0.4])
    xlabel('time (seconds)')
    ylabel('meters')
    legend()

    sgtitle('Traverse trajectories')

end