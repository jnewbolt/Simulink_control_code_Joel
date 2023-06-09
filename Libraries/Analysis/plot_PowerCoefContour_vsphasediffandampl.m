close all;

% Plot vs phase
f_star_sorted = f_star_commanded;
U_star_sorted = 1./f_star_sorted;
% Repeat the first row because data is periodic -180deg = 180deg
phase12_sorted = [phase12;-phase12(1,:)];
A_star_sorted = [A_star_measured;A_star_measured(1,:)];
powercoef_mean_sorted = [powercoef_mean;powercoef_mean(1,:)];

% Plot acceleration limit
acc_limit = 4.9; % Acceleration limit in m/s^2
v_limit = 0.5;

hold on
%% Plot Cp vs. A* and phase12
contourf(phase12_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none')
caxis([-1 0.2])
colormap(bluewhitered)

contour(phase12_sorted,A_star_sorted,powercoef_mean_sorted,[-1e-6 -1e-6],'LineWidth',4,'LineColor','k','LineStyle','--')
scatter(phase12_sorted,A_star_sorted,60,'.','k')
grid on
xlabel('Phase between foil and vibrissae (degrees)')
ylabel('{\it A} * = {\it A/d}')
xlim([-190 180])
ylim([-0.04 1.12])
set(gca, 'Layer', 'top')
grid off

c=colorbar();
c.Label.String = '{\it C}_P';
% c.Label.Interpreter = 'Latex';
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off