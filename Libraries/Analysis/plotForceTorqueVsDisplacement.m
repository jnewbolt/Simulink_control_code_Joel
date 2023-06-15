function plotForceTorqueVsDisplacement(time_star,T,freq,pitch_measured,heave_star_measured,liftcoef,...
        torquezcoef,titlePlots)
%(time_star,pitch_measured,heave_star_measured,liftcoef,dragcoef, ...
%    powercoef,torqueliftcoef,torquedragcoef,torquezcoef,num_cyc,plotTitle)
figure
hold on
grid on
plot(heave_star_measured,liftcoef,'DisplayName','{\it C}_L vs. {\it y/d}','Color','red','LineWidth',4)
plot(pitch_measured,torquezcoef,'DisplayName','{\it C}_{z\tau} vs. {\it \theta/\pi}','Color',"#FF00FF",'LineWidth',4)
quiver(heave_star_measured(1),liftcoef(1),heave_star_measured(2)-heave_star_measured(1),liftcoef(2)-liftcoef(1),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')
halfPeriodStep = round(1/(2*freq*T));
quiver(heave_star_measured(halfPeriodStep),liftcoef(halfPeriodStep),...
    heave_star_measured(halfPeriodStep+1)-heave_star_measured(halfPeriodStep),liftcoef(halfPeriodStep+1)-liftcoef(halfPeriodStep),...
    'LineWidth',4,'MaxHeadSize',1,'Color','black','AutoScaleFactor',50,'HandleVisibility','off')%
% plot(time_star,zeros(1,length(time_star)),'HandleVisibility','off','Color','#333333','LineWidth',3)
% plot(time_star,pitch_measured,'DisplayName','{\it \theta/\pi } ','Color','black','LineWidth',4,'LineStyle','--')
% plot(time_star,heave_star_measured,'DisplayName','{\it y/d}','Color','black','LineWidth',4)
% plot(time_star,liftcoef,'DisplayName','{\it C}_L','Color','red','LineWidth',4)
% plot(time_star,dragcoef,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% plot(time_star,torqueliftcoef,'DisplayName','{\it C}_{L\tau}','Color',"#A2142F",'LineWidth',4)
% plot(time_star,torquedragcoef,'DisplayName','{\it C}_{D\tau}','Color',"#0072BD",'LineWidth',4)
% plot(time_star,torquezcoef,'DisplayName','{\it C}_{z\tau}','Color',"#FF00FF",'LineWidth',4)
% plot(time_star,dragtorquecoef,'DisplayName','{\it C}_{\tau D}','Color','green','LineWidth',4)
% plot(time_star,powercoef,'DisplayName','C_P','Color','magenta')
hold off

legend()
ylim([-15 15])
xlim([-1 1])
xlabel('Displacement (dimensionless)')
ylabel('Force or torque coefficients')
title(titlePlots)
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
end
