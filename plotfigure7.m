getplotid
figure(7)



%% panel 2,4,6
if (VIR_ind ==4 && khyd_ind ==4)
    AMPK_AUC = AMPK_AUC./max(AMPK_AUC,[],'all');
    mTORC1_AUC = mTORC1_AUC./max(mTORC1_AUC,[],'all');
    mTORC2_AUC = mTORC2_AUC./max(mTORC2_AUC,[],'all');
subplot(3,2,2)
currentheatmap1 = heatmap(khyd_vary,parscan,AMPK_AUC);
ylabel('\bf{V_{IR} Multiplier}')
xlabel('\bf{Hydrolysis Rate Multiplier}')
title('Normalized AMPK AUC')
plotformatheatmap
currentheatmap1.FontSize = 12;
% currentheatmap1.ColorScaling ='log';
set(gca,'CellLabelColor','none')
caxis([0,1])

%% panel 2
subplot(3,2,4)
currentheatmap2 = heatmap(khyd_vary,parscan,mTORC1_AUC);
ylabel('\bf{V_{IR} Multiplier}')
xlabel('\bf{Hydrolysis Rate Multiplier}')
title('Normalized mTORC1 AUC')
plotformatheatmap
currentheatmap2.FontSize = 12;
% currentheatmap2.ColorScaling ='log';
set(gca,'CellLabelColor','none')
caxis([0,1])

%% panel 3
subplot(3,2,6)
currentheatmap3 = heatmap(khyd_vary,parscan,mTORC2_AUC);
ylabel('\bf{V_{IR} Multiplier}')
xlabel('\bf{Hydrolysis Rate Multiplier}')
title('Normalized mTORC2 AUC')
plotformatheatmap
caxis([0,1])

currentheatmap3.FontSize = 12;
% currentheatmap3.ColorScaling ='log';
set(gca,'CellLabelColor','none')
colormap(colorcet('L18'))

end
%% panel 1,3,5
subplot(3,2,1)
hold on
plot(t,y(:,pAMPK),'Color',colorvector(VIR_ind,:),'LineStyle',char(markervector(khyd_ind)))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
plotformat
xlim([-50,100])
xticklabels('auto')
title('AMPK')


subplot(3,2,3)
hold on
plot(t,y(:,pmTORC1),'Color',colorvector(VIR_ind,:),'LineStyle',char(markervector(khyd_ind)))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
plotformat
xlim([-50,100])
xticklabels('auto')
title('mTORC1')

subplot(3,2,5)
hold on
plot(t,y(:,pmTORC2),'Color',colorvector(VIR_ind,:),'LineStyle',char(markervector(khyd_ind)))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
plotformat
xlim([-50,100])
xticklabels('auto')
title('mTORC2')
