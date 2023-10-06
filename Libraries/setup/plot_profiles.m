%% plot profiles

function plot_profiles(P,arg1, arg2, arg3, arg4, arg5)

    if size(arg1,2) == 5
    % trajectories all in one array including reference signal
        data_p1 = arg1(:,1);
        data_h1 = arg1(:,2);
        data_p2 = arg1(:,3);
        data_h2 = arg1(:,4);
        data_ref = arg1(:,5);
    elseif size(arg1,2) == 4
    % trajectories all in one array excluding reference signal
        data_p1 = arg1(:,1);
        data_h1 = arg1(:,2);
        data_p2 = arg1(:,3);
        data_h2 = arg1(:,4);
        data_ref = 0;
    else
    % trajectories in differente individual vectors
        data_p1 = arg1;
        data_h1 = arg2;
        data_p2 = arg3;
        data_h2 = arg4;
        if ~exist('arg5', 'var')
            data_ref = 0;
        else
        data_ref = arg5;
        end
    end
    timeSeconds = linspace(1,length(data_p1),length(data_p1))'/P.sampleRate;
    
    figure
    subplot(2,1,1)
    plot(timeSeconds,data_p1); hold on;
    plot(timeSeconds,data_p2);
    plot(timeSeconds,data_ref*270,'r--'); hold off;
    ylim([-90,270])
    xlabel('time (seconds)')
    ylabel('degrees')
    legend('pitch 1', 'pitch 2', 'ref')

    subplot(2,1,2)
    plot(timeSeconds,data_h1); hold on;
    plot(timeSeconds,data_h2);
    plot(timeSeconds,data_ref,'r--'); hold off;
    ylim([-0.15,0.4])
    xlabel('time (seconds)')
    ylabel('meters')
    legend('heave 1', 'heave 2', 'ref')

    sgtitle('Traverse trajectories')

end