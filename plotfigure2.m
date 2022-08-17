getplotid
figure(22)


%% panel 1
subplot(5,3,1)
hold on
plot(t,y(:,Glut ))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Glutamate')
plotformat
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
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
xlim([-100,2000])
xticklabels('auto')
yticklabels('auto')

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



