%% Tandem foil analysis
% Force and efficiency calculation for tandem foil experiments
% 2023/04/17

clear;

addpath(genpath('Libraries'))

foil = foils_database('A3E');
fred = 0.12;

P1star_vec = [ 40    50    70 ]; lP1 = length(P1star_vec);
P2star_vec = 70; lP2 = length(P2star_vec);
H2star_vec = [ 0.4000    0.6000    0.8000    1.0000    1.2000 ]; lH2 = length(H2star_vec);
H2_vec = H2star_vec*foil.chord;
phase_vec = [ -180  -120   -60     0    60   120   180 ]; lph = length(phase_vec);

eff1 = nan(lP1,lP2,lH2,lph);
eff1_std = nan(lP1,lP2,lH2,lph);
eff2 = nan(lP1,lP2,lH2,lph);
sff2_std = nan(lP1,lP2,lH2,lph);
eff_sys = nan(lP1,lP2,lH2,lph);
eff_sys_std = nan(lP1,lP2,lH2,lph);

for ii = 1:lP1
    for jj = 1:lP2
        for kk = 1:lH2
            for qq = 1:lph
                aT4 = atan(-2*pi*(H2star_vec(kk))*fred) + deg2rad(P1star_vec(ii));
                FOLDERNAME = '\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20230414_TandemFriday_4c_separation_3alphaSweep_APHPH_A3E\';
                FILENAME = ['20230414_TandemThursday_4c_separation_3alphaSweep_',...
                    'aT4=',num2str(aT4),'_p2=',num2str(P1star_vec(ii)),'deg_h2=',num2str(H2_vec(kk)),'c_ph=',num2str(phase_vec(qq)),'deg.mat'];
                load(fullfile(FOLDERNAME,FILENAME));

                [kinematics, forces, nonDimForces, results] = calculateEfficiency(experiment, foil, out, freq, srate, 'filter', 20);

                eff1(ii,jj,kk,qq) = results.eff1;
                eff1_std(ii,jj,kk,qq) = results.eff1_std;
                eff2(ii,jj,kk,qq) = results.eff2;
                sff2_std(ii,jj,kk,qq) = results.eff2_std;
                eff_sys(ii,jj,kk,qq) = results.eff_sys;
                eff_sys_std(ii,jj,kk,qq) = results.eff_sys_std;

            end
        end
    end
end

savelocation = '\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20230414_TandemFriday_4c_3alphaSweep_efficiency_results';
finalfilename = ('20230417_TandemFriday_4c_separation_3alphaSweep_15_33_68_PHPh.mat');

save(fullfile(savelocation,finalfilename),'eff1','eff1_std','eff2','eff2_std','eff_sys','eff_sys_std');

