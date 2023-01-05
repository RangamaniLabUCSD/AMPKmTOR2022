getplotid
figure(3)




%% panel 2
subplot(3,2,1)
hold on
plot(t,y(:,Ca_C))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Cytosolic Calcium')
plotformat
xlim([-100,200])
xline(0,'-',{'1 Hz Glutamate Stimulus'},'LabelHorizontalAlignment' , {'left'})
xline(50,'-',{'Stimulus End'},'LabelHorizontalAlignment' , {'right'})

xticklabels('auto')
yticklabels('auto')

%% panel 1
subplot(3,2,2)
hold on
plot(t,y(:,AMP)./y(:,ATP))
xlabel('Time [s]')
ylabel('Concentration Ratio')
title('AMP/ATP Ratio')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 3
subplot(3,2,3)
hold on
plot(t,y(:,pAMPK))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pAMPK')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 4
subplot(3,2,4)
hold on
plot(t,y(:,pmTORC1))
%plot(t,y(:,pmTORC1)./y(:,mTORC1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC1')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 5
subplot(3,2,5)
hold on
plot(t,y(:,pmTORC2))
%plot(t,y(:,pmTORC2)./y(:,mTORC2))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC2')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 6
subplot(3,2,6)
hold on
plot(t,y(:,pULK1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pULK1')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

% barlabel = categorical({'Calcium','pAMPK','pmTORC1','pmTORC2'});
% bar(barlabel,[CaC_AUC./max(CaC_AUC);AMPK_AUC./max(AMPK_AUC);mTORC1_AUC./max(mTORC1_AUC);mTORC2_AUC./max(mTORC2_AUC)])
% ylabel('Normalized Area Under Curve')
% legend('0.1 Hz','1','10','100')
% 
% %% panel 7 AUC
% subplot(3,3,7)
%  barlabel = categorical({'pAMPK','pmTORC1','pmTORC2'});
%  bar(barlabel,[normAMPK(1,:);normmTOR1(1,:);normmTOR2(1,:)])
%  ylabel('Relative Area Under Curve')
%  legend('1 Hz')
% ylabel('Normalized AUC')
% ylim([0.9,1.2])
% plotformat
% 
% %% panel 8 Period
% subplot(3,3,8)
%  bar(barlabel,[period_a;period_m1;period_m2])
% ylim([15,40])
%  ylabel('Period [s]')
% 
% plotformat
% %% panel 9 maxamplitude
% subplot(3,3,9)
%  bar(barlabel,100.*[maxamp_a;maxamp_m1;maxamp_m2])
%  ylabel('Amplitude [% Change]')
% plotformat
% 
% % %% panel 4 ttp
% subplot(2,2,4)
%  bar(barlabel,[CaC_AUC;AMPK_AUC;mTORC1_AUC;mTORC2_AUC])
%  ylabel('Normalized Area Under Curve')
%  legend('0.1 Hz','1','10','100')
% xlabel('Species')
% ylabel('Normalized AUC')

plotformat
% %% panel 7
% subplot(3,3,7)
% hold on
% plot(t,y(:,))
% xlabel('Time [s]')
% ylabel('Concentration [\muM]')
% title('')
% plotformat
% 
% 
% 
% 
% 
% %% panel 8
% subplot(3,3,8)
% hold on
% plot(t,y(:,))
% xlabel('Time [s]')
% ylabel('Concentration [\muM]')
% title('')
% plotformat
% 
% %% panel 9
% subplot(3,3,9)
% hold on
% plot(t,y(:,))
% xlabel('Time [s]')
% ylabel('Concentration [\muM]')
% title('')
% plotformat
% 




getplotid
figure(23)


%% panel 1
subplot(5,3,1)
hold on
plot(t,y(:,Glut ))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Glutamate')
plotformat
xlim([-100,200])
xline(0,'-',{'1 Hz Glutamate Stimulus'},'LabelHorizontalAlignment' , {'left'})
xline(50,'-',{'Stimulus End'},'LabelHorizontalAlignment' , {'right'})
xticklabels('auto')
yticklabels('auto')

%% panel 2

subplot(5,3,2)
hold on
plot(t,y(:,Ca_C))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Cytosolic Calcium')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 3
subplot(5,3,3)
hold on
plot(t,y(:,NMDA_O))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Open NMDA')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 4
subplot(5,3,4)
hold on
plot(t,y(:,AMPA_O))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Open AMPA')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 5
subplot(5,3,5)
hold on
%plot(t,y(:,pmTORC1)./y(:,mTORC1))
plot(t,y(:,DIM))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Active mGLUR')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 6
subplot(5,3,6)
hold on
plot(t,y(:,IP3))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('IP3')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 7
subplot(5,3,7)
hold on
plot(t,y(:,CAMKK2))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('CAMKK2')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 8
subplot(5,3,8)
hold on
plot(t,y(:,AMP)./y(:,ATP))
xlabel('Time [s]')
ylabel('Concentration Ratio')
title('AMP/ATP')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')


%% panel 9
subplot(5,3,9)
hold on
plot(t,y(:,pAMPK))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pAMPK')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 10
subplot(5,3,10)
hold on
plot(t,y(:,pmTORC1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC1')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 11
subplot(5,3,11)
hold on
plot(t,y(:,pmTORC2))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC2')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 12
subplot(5,3,12)
hold on
plot(t,y(:,pAKT))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pAKT')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 13
subplot(5,3,13)
hold on
plot(t,y(:,pIRS))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Active IRS')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 14
subplot(5,3,14)
hold on
plot(t,y(:,pULK1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pULK1')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 15
subplot(5,3,15)
hold on
plot(t,y(:,SIRT1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('SIRT1')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')
