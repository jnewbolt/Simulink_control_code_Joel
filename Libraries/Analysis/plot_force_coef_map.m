function plot_force_coef_map(pitchAmpDegVec,heaveAmpMetersVec,forceCoefMeans,angleOfAttackMaxDegMatrix,plotTitle,colorbarLabel)

figure
hold on
% Plot Cp vs. A* and phase12
% powerCoefMeanReshape = reshape(powerCoefMeanG,length(phaseLagWbehindGvec),length(heaveAmpMetersGvec));
contourf(pitchAmpDegVec,heaveAmpMetersVec,forceCoefMeans,120,'LineStyle','none')
clim([-0 3])
colormap(bluewhitered)

contour(pitchAmpDegVec,heaveAmpMetersVec,forceCoefMeans,[-1e-6 -1e-6],'LineWidth',4,'LineColor','k','LineStyle','--')
% [C,h] = contour(pitchAmpDegVec,heaveAmpMetersVec,angleOfAttackMaxDegMatrix,[-20 -10 0 10 20],'LineWidth',3,'LineColor','k','ShowText','on');
% clabel(C,h,'FontSize',20)
% scatter(phase_vec,H2star_vec,60,'.','k')
grid on
xlabel('Foil pitch amplitude (degrees)')
ylabel('{\it A} * = {\it A/c}')
xlim([0 30])
ylim([0 0.5])
set(gca, 'Layer', 'top')
grid off

c=colorbar();
c.Label.String = colorbarLabel;
% c.Label.Interpreter = 'Latex';
% plotTitle = ['Object behind foil, foil pitch +/-',num2str(pitchAmpDegG,2),'deg and +/-',num2str(heaveAmpMetersG/chordMetersG,2),'chord'];
title(plotTitle);
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off
end

