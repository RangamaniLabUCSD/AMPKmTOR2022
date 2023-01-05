getplotid
figure(6)

% barlabel = categorical({'Calcium','pAMPK','pmTORC1','pmTORC2'});
% bar(barlabel,[CaC_AUC./max(CaC_AUC);AMPK_AUC./max(AMPK_AUC);mTORC1_AUC./max(mTORC1_AUC);mTORC2_AUC./max(mTORC2_AUC)])
% ylabel('Normalized Area Under Curve')
% legend('0.1 Hz','1','10','100')

%panel 6 SS
subplot(3,3,6)
 barlabel = categorical({'pAMPK','pmTORC1','pmTORC2'});
 bar(barlabel,[AMPK_SS(1,:);mTOR1_SS(1,:);mTOR2_SS(1,:)])
 ylabel('Normalized Steady State')

plotformat



%% panel 1 AUC
subplot(3,3,7)
 barlabel = categorical({'pAMPK','pmTORC1','pmTORC2'});
 bar(barlabel,[normAMPK;normmTOR1;normmTOR2])
 ylabel('Relative Area Under Curve')
 legend('0.1 \times k_{HYD}','1\times k_{HYD}','2\times k_{HYD}','5\times k_{HYD}')
ylabel('Normalized AUC')
plotformat

%% panel 2 Time to Equilibrium
subplot(3,3,8)
 bar(barlabel,[tte_a;tte_m1;tte_m2])
% ylim([0,200])
 ylabel('Time to Equilibrium [s]')

plotformat
%% panel 3 maxamplitude
subplot(3,3,9)
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
subplot(3,3,1)
 legend('0.1 \times k_{HYD}','1\times k_{HYD}','2\times k_{HYD}','5\times k_{HYD}')