function plotForceTorqueDisplacementVsTime(time_star,pitch_measured,heave_star_measured,liftcoef,dragcoef, ...
    powercoef,torqueliftcoef,torquedragcoef,torquezcoef,num_cyc,plotTitle)
figure('Position', [100 80 1600 900])
grid on
hold on

plot(time_star,zeros(1,length(time_star)),'HandleVisibility','off','Color','#333333','LineWidth',3) % emphasize y=0
% Plot displacements
plot(time_star,pitch_measured,'DisplayName','{\it \theta } (rad)','Color','black','LineWidth',4,'LineStyle','--')
plot(time_star,heave_star_measured,'DisplayName','{\it y/d}','Color','black','LineWidth',4)
% Plot force coefficients
plot(time_star,liftcoef,'DisplayName','{\it C}_L','Color','red','LineWidth',4)
plot(time_star,dragcoef,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
plot(time_star,torqueliftcoef,'DisplayName','{\it C}_{L\tau}','Color',"#A2142F",'LineWidth',4)
plot(time_star,torquedragcoef,'DisplayName','{\it C}_{D\tau}','Color',"#0072BD",'LineWidth',4)
plot(time_star,torquezcoef,'DisplayName','{\it C}_{z\tau}','Color',"#FF00FF",'LineWidth',4)
% Plot average values
plot(time_star, mean(liftcoef)*ones(length(time_star),1),'Color','red','LineWidth',4,'LineStyle','-.','DisplayName','mean {\it C}_L');
plot(time_star, mean(dragcoef)*ones(length(time_star),1),'Color','blue','LineWidth',4,'LineStyle','-.','DisplayName','mean {\it C}_D');

hold off

legend()
ylim([-5 5])
xlim([9.75 12.75])
xlabel('Time (cycles)')
ylabel('Displacements and force coefficients')
title(plotTitle)
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
end
