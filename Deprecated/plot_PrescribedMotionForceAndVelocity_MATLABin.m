function plotForceAndPosition(time_star,pitch_measured,heave_star_measured,liftcoef,dragcoef,powercoef,num_cyc,dragtorquecoef...
    ,plotTitle)
figure
hold on
grid on
plot(time_star,pitch_measured,'DisplayName','{\it \theta}','Color','black','LineWidth',4,'LineStyle','--')
plot(time_star,heave_star_measured,'DisplayName','{\it y/d}','Color','black','LineWidth',4)
plot(time_star,liftcoef,'DisplayName','{\it C}_L','Color','red','LineWidth',4)
plot(time_star,dragcoef,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% plot(time_star,dragtorquecoef,'DisplayName','{\it C}_{\tau D}','Color','green','LineWidth',4)
% plot(time_star,powercoef,'DisplayName','C_P','Color','magenta')
hold off

legend()
ylim([-15 15])
xlim([9.75 12.75])
xlabel('Time (cycles)')
title(plotTitle)
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
end
