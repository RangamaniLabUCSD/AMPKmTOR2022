%Plot 3D bifurcation diagrams for pmTORC1 and pmTORC2 varying V_IR
%and AMPK from saved matrices from simulation
%need to run bifurcation_2 first

%Import saved data for plotting
%inform user to run bifurcationAnalysis1_3D.m
data = open('bifVIR_Khyd_1.mat'); %Import saved simulation data
ypar = 'K_{HYD}';
%parval =data.kHYD;
Abundance =data.Abundance; %Define AMPK abundance
max_val = data.output_max_val; %Get matrix containing the maximum values
min_val = data.output_min_val; %Get matrix containing the minimum values
sum_val = data.output_sum; %Get matrix containing the sum
product_val = data.output_product; %Get matrix containing the product
div_val =  data.output_div; %Get matrix containing the quotient
dim = size(max_val); %Find the dimensions of the matrix

%Define parameter range for V_IR
param_range = data.param_range;
V_IR = data.param_range;
lim = [0 max(V_IR)]; %Define limit of V_IR to display
view_angle = [-30, 45]; %Define view angle of plot

%Define mesh for plotting scatter plots
[X, Y] = meshgrid(V_IR, Abundance);
X =  reshape(X, [dim(1)*dim(2), 1]);
Y =  reshape(Y, [dim(1)*dim(2), 1]);

%Reshape matrices
max_val_new = reshape(max_val, [dim(1)*dim(2), 2]);
min_val_new = reshape(min_val, [dim(1)*dim(2), 2]);
sum_val_new = reshape(sum_val, [dim(1)*dim(2), 2]);
product_val_new = reshape(product_val, [dim(1)*dim(2), 2]);
div_val_new = reshape(div_val, [dim(1)*dim(2), 2]);

yrange = [min(Abundance), max(Abundance) ] ;

%Plot 3D surface plots
%pmTORC1 
figure(1)
subplot(3,2,5)
surf(V_IR,Abundance,reshape(max_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none','edgecolor','none')
hold on 
surf(V_IR,Abundance,reshape(min_val(:, 1:length(V_IR), 1), [dim(1), dim(2)]), 'LineStyle','none','edgecolor','none')
scatter3(X, Y, max_val_new(:,1), 1, 'k.')
scatter3(X, Y, min_val_new(:,1), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel(ypar)
zlabel('pmTORC1')
xlim(lim)
view(view_angle)
plotformat

%pmTORC2
subplot(3,2,6)
surf(V_IR,Abundance,reshape(max_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
hold on 
surf(V_IR,Abundance,reshape(min_val(:, 1:length(V_IR), 2), [dim(1), dim(2)]), 'LineStyle','none')
scatter3(X, Y, max_val_new(:,2), 1, 'k.')
scatter3(X, Y, min_val_new(:,2), 1, 'k.')
hold off
xlabel('V\_IR')
ylabel(ypar)
zlabel('pmTORC2')
xlim(lim)
view(view_angle)
plotformat

xticklabels({0,5,10,15,20})
subplot(3,2,3)
plot1 = min_val(1,:,1);
plot2 = max_val(1,:,1);
plot(param_range,plot1,'b')
hold on
plot(param_range,plot2,'b')

plot1 = min_val(150,:,1);
plot2 = max_val(150,:,1);
plot(param_range,plot1,'r')
hold on
plot(param_range,plot2,'r')
ylabel('Steady State [\muM]')
xlabel('V\_IR')
title('mTORC1')

plotformat

subplot(3,2,4)
plot1 = min_val(1,:,2);
plot2 = max_val(1,:,2);
plot(param_range,plot1,'b')
hold on
plot(param_range,plot2,'b')
ylabel('Steady State [\muM]')
xlabel('V\_IR')
title('mTORC2')
plotformat


plot1 = min_val(150,:,2);
plot2 = max_val(150,:,2);
plot(param_range,plot1,'r')
hold on
plot(param_range,plot2,'r')
plotformat
legend(num2str(Abundance(1)),'',num2str(Abundance(150)),'');
leg=legend('show');


subplot(3,2,1)
hold on
plot1 = min_val(:,30,1);
plot2 = max_val(:,30,1);
plot(Abundance,plot1,'Color',[0.9290 0.6940 0.1250])
hold on
plot(Abundance,plot2,'Color',[0.9290 0.6940 0.1250])
ylabel('Steady State [\muM]')
xlabel('k_{hyd}')
title('mTORC2')
plotformat

plot1 = min_val(:,1,1);
plot2 = max_val(:,1,1);
plot(Abundance,plot1,'Color',[0.4940 0.1840 0.5560])
hold on
plot(Abundance,plot2, 'Color',[0.4940 0.1840 0.5560])
ylabel('Steady State [\muM]')
xlabel('k_{hyd}')
title('mTORC2')
plotformat
legend(num2str(param_range(30)),'',num2str(param_range(1)),'');
leg=legend('show');


subplot(3,2,2)


hold on
plot1 = min_val(:,30,1);
plot2 = max_val(:,30,1);
plot(Abundance,plot1,'Color',[0.9290 0.6940 0.1250])
hold on
plot(Abundance,plot2,'Color',[0.9290 0.6940 0.1250])
ylabel('Steady State [\muM]')
xlabel('k_{hyd}')
title('mTORC1')
plotformat

plot1 = min_val(:,1,2);
plot2 = max_val(:,1,2);
plot(Abundance,plot1,'Color',[0.4940 0.1840 0.5560])
hold on
plot(Abundance,plot2, 'Color',[0.4940 0.1840 0.5560])
ylabel('Steady State [\muM]')
xlabel('k_{hyd}')
title('mTORC2')
plotformat
legend(num2str(param_range(30)),'',num2str(param_range(1)),'');
leg=legend('show');