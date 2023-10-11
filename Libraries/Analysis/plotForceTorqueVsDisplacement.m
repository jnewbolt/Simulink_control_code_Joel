function plotForceTorqueVsDisplacement(timeStar,sampleRate,freq,pitchMeasured,heaveStarMeasured,liftCoef,...
        torquezCoef,titlePlots)

figure('Position', [100 80 1600 900])
hold on
grid on
plot(heaveStarMeasured,liftCoef,'DisplayName','{\it C}_L vs. {\it y/d}','Color','red','LineWidth',4)
plot(pitchMeasured,torquezCoef,'DisplayName','{\it C}_{z\tau} vs. {\it \theta} (rad)','Color',"#FF00FF",'LineWidth',4)
% Show arrows to indicate the direction of the force trajectory around the force loop with time
quiverIndex= round(max(timeStar/2))*sampleRate/freq;
quiver(heaveStarMeasured(quiverIndex),liftCoef(quiverIndex),heaveStarMeasured(quiverIndex+1)-heaveStarMeasured(quiverIndex),liftCoef(quiverIndex+1)-liftCoef(quiverIndex),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')
% annotation('arrow',[heaveStarMeasured(quiverIndex) heaveStarMeasured(quiverIndex+1)]/1.2+0.5,[liftCoef(quiverIndex) liftCoef(quiverIndex+1)]/16+0.5)
halfPeriodStep = quiverIndex+round(sampleRate/(2*freq));
quiver(heaveStarMeasured(halfPeriodStep),liftCoef(halfPeriodStep),...
    heaveStarMeasured(halfPeriodStep+1)-heaveStarMeasured(halfPeriodStep),liftCoef(halfPeriodStep+1)-liftCoef(halfPeriodStep),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')%
% Show arrows to indicate the direction of the torque trajectory around the force loop with time
quiver(pitchMeasured(quiverIndex),torquezCoef(quiverIndex),pitchMeasured(quiverIndex+1)-pitchMeasured(quiverIndex),torquezCoef(quiverIndex+1)-torquezCoef(quiverIndex),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')
halfPeriodStep = round(sampleRate/(2*freq));
quiver(pitchMeasured(halfPeriodStep),torquezCoef(halfPeriodStep),...
    pitchMeasured(halfPeriodStep+1)-pitchMeasured(halfPeriodStep),torquezCoef(halfPeriodStep+1)-torquezCoef(halfPeriodStep),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')%
hold off

legend()
ylim([-8 8])
xlim([-0.6 0.6])
xlabel('Displacement (dimensionless)')
ylabel('Force or torque coefficients')
title(titlePlots)
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
end
