function plotForceTorqueVsDisplacement(timeStar,sampleRate,freq,pitch_measured,heave_star_measured,liftcoef,...
        torquezcoef,titlePlots)
%(time_star,pitch_measured,heave_star_measured,liftcoef,dragcoef, ...
%    powercoef,torqueliftcoef,torquedragcoef,torquezcoef,num_cyc,plotTitle)
figure('Position', [100 80 1600 900])
hold on
grid on
plot(heave_star_measured,liftcoef,'DisplayName','{\it C}_L vs. {\it y/d}','Color','red','LineWidth',4)
plot(pitch_measured,torquezcoef,'DisplayName','{\it C}_{z\tau} vs. {\it \theta} (rad)','Color',"#FF00FF",'LineWidth',4)
% Show arrows to indicate the direction of the force trajectory around the force loop with time
quiver(heave_star_measured(1),liftcoef(1),heave_star_measured(2)-heave_star_measured(1),liftcoef(2)-liftcoef(1),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')
halfPeriodStep = round(sampleRate/(2*freq));
quiver(heave_star_measured(halfPeriodStep),liftcoef(halfPeriodStep),...
    heave_star_measured(halfPeriodStep+1)-heave_star_measured(halfPeriodStep),liftcoef(halfPeriodStep+1)-liftcoef(halfPeriodStep),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')%
% Show arrows to indicate the direction of the torque trajectory around the force loop with time
quiver(pitch_measured(1),torquezcoef(1),pitch_measured(2)-pitch_measured(1),torquezcoef(2)-torquezcoef(1),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')
halfPeriodStep = round(sampleRate/(2*freq));
quiver(pitch_measured(halfPeriodStep),torquezcoef(halfPeriodStep),...
    pitch_measured(halfPeriodStep+1)-pitch_measured(halfPeriodStep),torquezcoef(halfPeriodStep+1)-torquezcoef(halfPeriodStep),...
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
