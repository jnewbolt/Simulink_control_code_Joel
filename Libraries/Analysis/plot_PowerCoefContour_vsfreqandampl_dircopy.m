close all;


f_star_grdpts = (0.1:0.02:0.26);
A_star_grdpts = (0.0:0.1:1.1);
[f_star_sorted,A_star_sorted] = meshgrid(f_star_grdpts,A_star_grdpts);
powercoef_mea
powercoef_mean_sorted = reshape(powercoef_mean,length(A_star_grdpts),length(f_star_grdpts));
% Sort the trials by frequency (needed for countorf function)
% [f_star_sorted,sort_index] = sortrows(f_star_commanded);
% U_star_sorted = 1./f_star_sorted;
% phase12_sorted = phase12(sort_index,:);
% A_star_sorted = A_star_measured(sort_index,:);
% powercoef_mean_sorted = powercoef_mean(sort_index,:);
% delay_sorted = delay(sort_index,:);
% powercoef_conv = powercoef_convtest(1,sort_index,:);

% Plot acceleration limit
% acc_limit = 6.5; % Acceleration limit in m/s^2
% v_limit = 0.3;

hold on

% Plot Cp vs. A* and f*
contourf(f_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% contourf(f_star_sorted,A_star_sorted,squeeze(powercoef_conv(1,:,:))./powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% quiver(f_star_sorted,A_star_sorted,gradA,gradf)

caxis([-1 0.2])
% colorbarpwn(-6.0,0.2,'colorN',[0 0.5 1],'log',1.5)
colormap(bluewhitered)

contour(f_star_sorted,A_star_sorted,powercoef_mean_sorted,[-1e-6 -1e-6],'LineWidth',4,'LineColor','k','LineStyle','--')
scatter(f_star_sorted,A_star_sorted,60,'.','k')
grid on
xlabel('{\it f} * = {\it f D/U}')
ylabel('{\it A} * = {\it A/D}')
xlim([0.27 0.41])% xlim([0.09 0.27]) %xlim([0.047 0.145]) %
ylim([-0.04 1.12])
% xticks([0.10 0.14 0.18 0.22 0.26])
set(gca, 'Layer', 'top')
grid off
% fstarvalues_extended = (0:0.025:0.3);
a_limit_curve = heaveaccelcommandlimit./(thcknss*(2*pi*(U.*fstarvector/thcknss)).^2);
v_limit_curve = heavevelocommandlimit./(thcknss*(2*pi*(U.*fstarvector/thcknss)));
% plot(fstarvector,a_limit_curve,'LineWidth',4,'Color','b','LineStyle','-.')
plot(fstarvector,v_limit_curve,'LineWidth',4,'Color','b','LineStyle','-.')

% % Plot Cp vs. A* and U*
% % contourf(U_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% % contourf(f_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% % quiver(f_star_sorted,A_star_sorted,gradA,gradf)
% caxis([-1.0 0.2])
% colormap(bluewhitered)
% contour(U_star_sorted,A_star_sorted,powercoef_mean_sorted,[0,0],'LineWidth',3,'LineColor','k','LineStyle','--')
% % scatter(U_star_sorted,A_star_sorted,[],'k')
% grid on
% xlabel('{\it U}* = {\it U/Df }')
% ylabel('{\it A}* = {\it A/D}')
% xlim([3 20])
% ylim([0 1.2])
% Ustarvalues_extended = (0:1:11);
% a_limit_curve = acc_limit./(chord*(2*pi*(flowspeed_fixed./(chord*Ustarvalues_extended))).^2);
% v_limit_curve = v_limit./(chord*(2*pi*(flowspeed_fixed./(chord*Ustarvalues_extended))));
% plot(Ustarvalues_extended,a_limit_curve)
% plot(Ustarvalues_extended,v_limit_curve)

c=colorbar();
c.Label.String = '{\it C}_P';
% c.Label.Interpreter = 'Latex';
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off