if exist('subplotvar','var') ==0  
subplotvar =1;
end

getplotid
figure(4)

%% panel 1
subplot(3,4,5)
hold on
plot(t,y(:,Ca_C))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('Cytosolic Calcium')
plotformat
xlim([-10,60])
xline(0,'-',{'Glutamate Stimulus'},'LabelHorizontalAlignment' , {'left'})
xline(5,'-',{'Stimulus End'},'LabelHorizontalAlignment' , {'right'})
xticklabels('auto')
yticklabels('auto')


%% panel 2
subplot(3,4,6)
hold on
plot(t,y(:,AMP)./y(:,ATP))
xlabel('Time [s]')
ylabel('Concentration Ratio')
title('AMP/ATP Ratio')
plotformat
xlim([-10,60])
xticklabels('auto')
yticklabels('auto')



%% panel 3
subplot(3,4,7)
hold on
plot(t,y(:,pAMPK))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pAMPK')
plotformat
xlim([-10,60])

xticklabels('auto')
yticklabels('auto')

%% panel 4
subplot(3,4,9)
hold on
plot(t,y(:,pmTORC1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC1')
plotformat
xlim([-100,200])
xticklabels('auto')
yticklabels('auto')

%% panel 5
subplot(3,4,10)
hold on
plot(t,y(:,pmTORC2))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pmTORC2')
plotformat
xlim([-100,300])
xticklabels('auto')
yticklabels('auto')

%% panel 6
subplot(3,4,11)
hold on
plot(t,y(:,pULK1))
xlabel('Time [s]')
ylabel('Concentration [\muM]')
title('pULK1')
plotformat
xlim([-100,300])
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
figure(4)


if subplotvar ==1
subplot(3,4,1)    
hold on
xlabel('Time [s]')
ylabel('Concentration [\muM]')
h=plot(t,y(:,12));
title('0.1 Hz and 1 Hz')
set(h, 'Color', '#0072BD')
end

if subplotvar ==2
    subplot(3,4,1)   
hold on
xlabel('Time [s]')
ylabel('Concentration [\muM]')
h=plot(t,y(:,12),'--');
title('0.1 Hz and 1 Hz')
set(h, 'Color', '#D95319')
legend('0.1', '1')
end

if subplotvar ==3
    subplot(3,4,2)   
hold on
xlabel('Time [s]')
ylabel('Concentration [\muM]')
h=plot(t,y(:,12));
title('10 Hz')
set(h, 'Color','#EDB120')
end

if subplotvar ==4
    subplot(3,4,3)   
hold on
xlabel('Time [s]')
ylabel('Concentration [\muM]')
h=plot(t,y(:,12));
title('100 Hz')
set(h, 'Color', '#7E2F8E')

end


plotformat
xlim([-05,10])
xticklabels('auto')
yticklabels('auto')
subplotvar = subplotvar+1;
if subplotvar>4
    subplotvar = 1;
end

