getplotid
figure(6)


%% panel 1
subplot(3,3,1)
hold on
plot(t,y(:,AMP)./y(:,ATP))
xlabel('Time [s]')
ylabel('Concentration Ratio')
title('AMP/ATP Ratio')
plotformat
xlim([-50,60])
xline(0,'-',{'Glutamate Stimulus'},'LabelHorizontalAlignment' , {'left'})
xline(5,'-',{'End'},'LabelHorizontalAlignment' , {'right'})
xticklabels('auto')
yticklabels('auto')
% 
% %% panel 2
% subplot(3,3,1)
% hold on
% plot(t,y(:,Ca_C))
% xlabel('Time [s]')
% ylabel('Concentration [\muM]')
% title('Cytosolic Calcium')
% plotformat
% xticklabels('auto')
% yticklabels('auto')
% xlim([-50,60])
% xline(0,'-',{'Glutamate Stimulus'},'LabelHorizontalAlignment' , {'left'})
% xline(5,'-',{'Stimulus End'},'LabelHorizontalAlignment' , {'right'})


%% panel 3
subplot(3,3,2)
hold on
plot(t,y(:,pAMPK))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pAMPK')
plotformat
xlim([-50,200])

xticklabels('auto')
yticklabels('auto')
%% panel 4
subplot(3,3,3)
hold on
plot(t,y(:,pmTORC1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC1')
xlim([-50,200])
plotformat
xticklabels('auto')
yticklabels('auto')
%% panel 5
subplot(3,3,4)
hold on
plot(t,y(:,pmTORC2))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC2')
xlim([-50,200])
plotformat
xticklabels('auto')
yticklabels('auto')
%% panel 6
subplot(3,3,5)
hold on
plot(t,y(:,pULK1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pULK1')
xlim([-50,200])
plotformat
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
% xticklabels('auto')
%yticklabels('auto')
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
% xticklabels('auto')
%yticklabels('auto')
% 
% 
% %% panel 9
% subplot(3,3,9)
% hold on
% plot(t,y(:,))
% xlabel('Time [s]')
% ylabel('Concentration [\muM]')
% title('')
% plotformat
% xticklabels('auto')
%yticklabels('auto')
% 
% 


