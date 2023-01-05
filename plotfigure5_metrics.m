getplotid
figure(4)

% barlabel = categorical({'Calcium','pAMPK','pmTORC1','pmTORC2'});
% bar(barlabel,[CaC_AUC./max(CaC_AUC);AMPK_AUC./max(AMPK_AUC);mTORC1_AUC./max(mTORC1_AUC);mTORC2_AUC./max(mTORC2_AUC)])
% ylabel('Normalized Area Under Curve')
% legend('0.1 Hz','1','10','100')

%% panel 1 AUC
subplot(3,4,4)
 barlabel = categorical({'pAMPK','pmTORC1','pmTORC2'});
 bar(barlabel,[normAMPK(1,:);normmTOR1(1,:);normmTOR2(1,:)])
 ylabel('Relative Area Under Curve')
 legend('0.1 Hz','1','10','100')
ylabel('Normalized AUC')
ylim([0.9,1.2])
plotformat

%% panel 2 Period
subplot(3,4,8)
 bar(barlabel,[tte_a;tte_m1;tte_m2])

 ylabel('Time to Equilibrium [s]')

plotformat
%% panel 3 maxamplitude
subplot(3,4,12)
 bar(barlabel,100.*[maxamp_a;maxamp_m1;maxamp_m2])
 ylabel('Amplitude [% Change]')
plotformat

% %% panel 4 ttp
% subplot(2,2,4)
%  bar(barlabel,[CaC_AUC;AMPK_AUC;mTORC1_AUC;mTORC2_AUC])
%  ylabel('Normalized Area Under Curve')
%  legend('0.1 Hz','1','10','100')
% xlabel('Species')
% ylabel('Normalized AUC')

plotformat