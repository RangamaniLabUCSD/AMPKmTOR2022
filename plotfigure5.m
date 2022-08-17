getplotid
figure(5)



%% panel 1
subplot(2,2,1)
currentheatmap1 = heatmap(glutscan,freq_scan,CaC_AUC./max(CaC_AUC));
xlabel('Glutamate Concentration [\muM]')
ylabel('Stimulus Frquency [Hz]')
title('Calcium AUC')
plotformatheatmap
currentheatmap1.FontSize = 12;
currentheatmap1.ColorScaling ='log';
set(gca,'CellLabelColor','none')

%% panel 2
subplot(2,2,2)
currentheatmap2 = heatmap(glutscan,freq_scan,AMPK_AUC./max(AMPK_AUC));
xlabel('Glutamate Concentration [\muM]')
ylabel('Stimulus Frquency [Hz]')
title('AMPK AUC')
plotformatheatmap
currentheatmap2.FontSize = 12;
currentheatmap2.ColorScaling ='log';
set(gca,'CellLabelColor','none')

%% panel 3
subplot(2,2,3)
currentheatmap3 = heatmap(glutscan,freq_scan,mTORC1_AUC./max(mTORC1_AUC));
xlabel('Glutamate Concentration [\muM]')
ylabel('Stimulus Frquency [Hz]')
title('mTORC1 AUC')
plotformatheatmap

currentheatmap3.FontSize = 12;
currentheatmap3.ColorScaling ='log';
set(gca,'CellLabelColor','none')

%% panel 4
subplot(2,2,4)
currentheatmap4 = heatmap(glutscan,freq_scan,mTORC2_AUC./max(mTORC2_AUC));
xlabel('Glutamate Concentration [\muM]')
ylabel('Stimulus Frquency [Hz]')
title('mTORC2 AUC')
plotformatheatmap

currentheatmap4.FontSize = 12;
currentheatmap4.ColorScaling ='log';
set(gca,'CellLabelColor','none')

