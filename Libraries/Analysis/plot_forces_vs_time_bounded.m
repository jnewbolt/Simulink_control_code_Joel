function plotForceTorqueDisplacementVsTime(timeStarPhaseAvg,pitchRadsPhaseAvg,pitchRadsPhaseStddev,heaveStarPhaseAvg,...
    heaveStarPhaseStddev,liftCoefPhaseAvg,liftCoefPhaseStddev,dragCoefPhaseAvg,dragCoefStddev,...
    powercoef,num_cyc,plotTitle)
figure('Position', [100 80 1600 900])
grid on
hold on

timeStarPhaseFill = [timeStarPhaseAvg'; flipud(timeStarPhaseAvg')];
pitchRadsPhaseAvgFill = [pitchRadsPhaseAvg+pitchRadsPhaseStddev;  flipud(pitchRadsPhaseAvg-pitchRadsPhaseStddev)];
heaveStarPhaseAvgFill = [heaveStarPhaseAvg+heaveStarPhaseStddev;  flipud(heaveStarPhaseAvg-heaveStarPhaseStddev)];
liftCoefPhaseAvgFill = [liftCoefPhaseAvg+liftCoefPhaseStddev;  flipud(liftCoefPhaseAvg-liftCoefPhaseStddev)];
dragCoefPhaseAvgFill = [dragCoefPhaseAvg+dragCoefStddev;  flipud(dragCoefPhaseAvg-dragCoefStddev)];
% plot(time_star,zeros(1,length(time_star)),'HandleVisibility','off','Color','#222222','LineWidth',3) % emphasize y=0
% Plot displacements
fill(timeStarPhaseFill,pitchRadsPhaseAvgFill,[0.6 0.6 0.6],'LineStyle','none','HandleVisibility','off')
plot(timeStarPhaseAvg,pitchRadsPhaseAvg,'DisplayName','{\it \theta } (rads)','Color',[0.3 0.3 0.3],'LineWidth',4)
fill(timeStarPhaseFill,heaveStarPhaseAvgFill,[0.3 0.3 0.3],'LineStyle','none','HandleVisibility','off')
plot(timeStarPhaseAvg,heaveStarPhaseAvg,'DisplayName','{\it y/c}','Color','black','LineWidth',4)
% Plot force coefficients
fill(timeStarPhaseFill,liftCoefPhaseAvgFill,[0.7 0 0],'LineStyle','none','HandleVisibility','off')
plot(timeStarPhaseAvg,liftCoefPhaseAvg,'DisplayName','{\it C}_L','Color','red','LineWidth',4)
fill(timeStarPhaseFill,dragCoefPhaseAvgFill,[0 0 0.7],'LineStyle','none','HandleVisibility','off')
plot(timeStarPhaseAvg,dragCoefPhaseAvg,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% plot(time_star,torqueliftcoef,'DisplayName','{\it C}_{L\tau}','Color',"#A2142F",'LineWidth',4)
% plot(time_star,torquedragcoef,'DisplayName','{\it C}_{D\tau}','Color',"#0072BD",'LineWidth',4)
% plot(time_star,torquezcoef,'DisplayName','{\it C}_{z\tau}','Color',"#FF00FF",'LineWidth',4)
% Plot average values
% plot(time_star, mean(liftcoef)*ones(length(time_star),1),'Color','red','LineWidth',4,'LineStyle','-.','DisplayName','mean {\it C}_L');
% plot(time_star, mean(dragcoef)*ones(length(time_star),1),'Color','blue','LineWidth',4,'LineStyle','-.','DisplayName','mean {\it C}_D');

hold off

lgd = legend();
lgd.Location = "southeast";
ylim([-3 3])
xlim([0 1])
xlabel('Time (cycles)')
ylabel({'Displacements and',' force coefficients'})
title(plotTitle)
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
end
