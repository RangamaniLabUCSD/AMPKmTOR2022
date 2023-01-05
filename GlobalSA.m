%Global Sensitivity Analysis for the mTORC1/AMPK/ULK1 Pathway
clc; clear;            
tic %start timing

%Set number of random samples
N = 10000; % Number of parameter sets

%Parameter ranges for Km, Vmax, kc and kd
Km_range = [1 1000];
Vmax_range = [1e-3 10] * 3600; %Convert units from per sec -> per h
Kc_range = [1e-6 1] * 3600; %Convert units from per sec -> per h
Kd_range = [1e-6 10] * 3600; %Convert units from per sec -> per h

tt = 0:500; %Timespan
init = leung_icsource ;
 load('finalpars.mat');
pulsetime = 0;
options = odeset('RelTol',1e-5, 'AbsTol',1e-8,  ...
    'InitialStep',1e-5, 'MaxOrder',5, 'BDF','on');
pars = finalpars;

 kHYD = pars(4);

metab_pars = [1:4,7:28];
cal_pars   = [29:60];
i_pars     = [103:163];


%Define parameters names for constructing the heatmap
parameter_names = {'KADP';'VmaxOP';'nH';'kHYD';'kstim';'kpost';'vAK';'kmt';'kmd';'kmm';'keqadk';'vmax20';...
'vmax21';'km20';'km21';'k12f';'k12r';'k13f';'k13r';'kaicar';'kAct';'kcatAMPK';...
'VforCK';'Kb';'Kia';'Kib';'Kiq';'KpCK';'KeqCK';'TCr';'k1';'b';...
'Ki';'Ka';'IRa';'k_-';'k';'VDAG';'kb';'ku';'Vp';'Kmp';...
'Vpkc';'Kmpkc';'kapkc';'kdpkc';'kplcI';'kplcD';'KmDAG';'VRYR';'Kar';'Kbr';...
'Kcr';'Kdr';'VSERCA';'Kp';'Rb';'Ru';'Rd';'Rr';'Ro';'Rc';...
'tauess';'tauesf';'taubss';'taubsf';'tdelay';'s';'V_rev';'N_NMDA';'BPAP max';'G_NMDA';...
'k_ext';'Kmext';'kdeg';'kpmleak';'k_buffcyt';'k_buffcyt_r';'k_ERBUFF';'k_ERBUFF_r';'kerleak';'G_AMPA';...
'kAMPA_1f';'kAMPA_2f';'kAMPA_3f';'kAMPA_4f';'kAMPA_5f';'kAMPA_6f';'kAMPA_7f';'kAMPA_8f';'kAMPA_1r';'kAMPA_2r';...
'kAMPA_3r';'kAMPA_4r';'kAMPA_5r';'kAMPA_6r';'kAMPA_7r';'kAMPA_8r';'k_cam_f';'k_cam_r';'Vm_kk2';'Km_kk2';...
'Ka_kk2';'kf_CaMDisc';'V_IR';'Km_IR';'V_pIR';'Km_pIR';'K_IRS_by_pIR';'Km_IRS_by_pIR';'V_pIRS';'Km_pIRS';...
'K_AKT_by_pIRS';'Km_AKT_by_pIRS';'K_AKT_by_pmTORC2';'Km_AKT_by_pmTORC2';'V_pAKT';'Km_pAKT';'K_mTORC1_by_pAKT';'Km_mTORC1_by_pAKT';'K_pmTORC1';'K_pmTORC1_by_pAMPK';...
'Km_pmTORC1_by_pAMPK';'K_pmTORC1_by_pULK1';'Km_pmTORC1_by_pULK1';'K_mTORC2_by_pIRS';'Km_mTORC2_by_pIRS';'K_mTORC2_by_pAMPK';'Km_mTORC2_by_pAMPK';'V_pmTORC2';'Km_pmTORC2';'K_DEPTOR_by_pmTORC1';...
'Km_DEPTOR_by_pmTORC1';'K_DEPTOR_by_pmTORC2';'Km_DEPTOR_by_pmTORC2';'V_pDEPTOR';'Km_pDEPTOR';'K_mTORC1_DEPTOR_form';'K_mTORC1_DEPTOR_diss';'K_mTORC2_DEPTOR_form';'K_mTORC2_DEPTOR_diss';'K_IRS_to_iIRS';...
'Km_IRS_to_iIRS';'V_iIRS';'Km_iIRS';'K_AMPK';'K_AMPK_by_SIRT1';'Km_AMPK';'K_pAMPK';'K_pAMPK_by_pULK1';'K_pAMPK_by_pmTORC1';'Km_pAMPK';...
'K_SIRT1';'K_SIRT1_by_pAMPK';'Km_SIRT1';'K_SIRT1_diss';'K_ULK1';'K_ULK1_by_pAMPK';'Km_ULK1';'K_pULK1';'K_pULK1_by_pmTORC1';'Km_pULK1';'AMPK_AMPact'}; 
        
parameter_names = strrep(parameter_names, '_', '\_'); % Replace '_' by '\_'
            
%Define variable names
variable_names = {'pAMPK','pAKT','pmTORC1','pmTORC2'};  

%Construct matrix for sampled parameter sets
%Each row of the array param corresponds to a sampled set of parameters
% param_input = {KADP;VmaxOP;nH;krest;kstim;kpost;vAK;kmt;kmd;kmm;keqadk;vmax20;...
% vmax21;km20;km21;k12f;k12r;k13f;k13r;kaicar;kAct;kcatAMPK;...
% VforCK;Kb;Kia;Kib;Kiq;KpCK;KeqCK;TCr;k1;b;...
% Ki;Ka;IRa;k_;k;VDAG;kb;ku;Vp;Kmp;...
% Vpkc;Kmpkc;kapkc;kdpkc;kplcI;kplcD;KmDAG;VRYR;Kar;Kbr;...
% Kcr;Kdr;VSERCA;Kp;Rb;Ru;Rd;Rr;Ro;Rc;...
% tauess;tauesf;taubss;taubsf;tdelay;s;V_rev;N_NMDA;BPAP_max;G_NMDA;...
% k_ext;Kmext;kdeg;kpmleak;k_buffcyt;k_buffcyt_r;k_ERBUFF;k_ERBUFF_r;kerleak;G_AMPA;...
% kAMPA_1f;kAMPA_2f;kAMPA_3f;kAMPA_4f;kAMPA_5f;kAMPA_6f;kAMPA_7f;kAMPA_8f;kAMPA_1r;kAMPA_2r;...
% kAMPA_3r;kAMPA_4r;kAMPA_5r;kAMPA_6r;kAMPA_7r;kAMPA_8r;k_cam_f;k_cam_r;Vm_kk2;Km_kk2;...
% Ka_kk2;kf_CaMDisc;V_IR;Km_IR;V_pIR;Km_pIR;K_IRS_by_pIR;Km_IRS_by_pIR;V_pIRS;Km_pIRS;...
% K_AKT_by_pIRS;Km_AKT_by_pIRS;K_AKT_by_pmTORC2;Km_AKT_by_pmTORC2;V_pAKT;Km_pAKT;K_mTORC1_by_pAKT;Km_mTORC1_by_pAKT;K_pmTORC1;K_pmTORC1_by_pAMPK;...
% Km_pmTORC1_by_pAMPK;K_pmTORC1_by_pULK1;Km_pmTORC1_by_pULK1;K_mTORC2_by_pIRS;Km_mTORC2_by_pIRS;K_mTORC2_by_pAMPK;Km_mTORC2_by_pAMPK;V_pmTORC2;Km_pmTORC2;K_DEPTOR_by_pmTORC1;...
% Km_DEPTOR_by_pmTORC1;K_DEPTOR_by_pmTORC2;Km_DEPTOR_by_pmTORC2;V_pDEPTOR;Km_pDEPTOR;K_mTORC1_DEPTOR_form;K_mTORC1_DEPTOR_diss;K_mTORC2_DEPTOR_form;K_mTORC2_DEPTOR_diss;K_IRS_to_iIRS;...
% Km_IRS_to_iIRS;V_iIRS;Km_iIRS;K_AMPK;K_AMPK_by_SIRT1;Km_AMPK;K_pAMPK;K_pAMPK_by_pULK1;K_pAMPK_by_pmTORC1;Km_pAMPK;...
% K_SIRT1;K_SIRT1_by_pAMPK;Km_SIRT1;K_SIRT1_diss;K_ULK1;K_ULK1_by_pAMPK;Km_ULK1;K_pULK1;K_pULK1_by_pmTORC1;Km_pULK1;}; 
variance_range = 0.2; % variance for the parameter values 0.5 = 50%
pars_low = pars.* (1-variance_range);
pars_high = pars.*(1 + variance_range);
%(Vmax_range(end)-Vmax_range(1))*lhsdesign(N,1) + Vmax_range(1)
for pari = 1:numel(pars)
param_input(pari,:)  =(pars_high(pari) - pars_low(pari)) * lhsdesign(N,1) + pars_low(pari);
end
param_input = param_input';

%Find number of parameters/variables
len_param = length(parameter_names); %Number of parameters
len_var = length(variable_names); %Number of variables

%Solve the system of ODEs for each sampled parameter set
tt = 0:5000; %Timespan

results = zeros(N, len_var); %Pre-allocate matrix to store results

%Solve the system of ODEs and store the terminal concentration vector in time 
parfor i = 1:N   
    kHYD = param_input(i,4);
    [t, temptraj] = ode15s( @leung2022_SERCA,tt,init, options,param_input(i,:), kHYD,pulsetime);
    results(i,:) = temptraj(end,[7,44,46,48]); %Capture terminal point of the simulation
    fprintf('i = %d\n',i) %Print current iteration number
end
fprintf('ODE solutions complete\n')

%Calculate mean and variance
results_mean = mean(results)';
results_var = var(results)';
results_table = table( variable_names', results_mean, results_var ); %Show table of mean / variance


%Plot the Probability Density Distributions for each Output
fig = figure(1);
hold on
for i = 1:len_var
   subplot(2,2,i)
   histogram(results(:,i), 250, 'Normalization', 'pdf');
   title(variable_names{i})
end
hold off

%Set plot properties
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Probability');
xlabel(han,'Abundance');
sgtitle('Probability Density Functions of the Output', 'fontsize', 10)


%Partial Rank Correlation Coefficient
%Find the PRCC for each output 
rho_final = zeros(len_param, len_var); %Pre-allocate matrix 
pvalues_final = zeros(len_param, len_var); %Pre-allocate matrix

%Compute the PRCC for each variable with respect to the model parameters
fprintf('PRCC solutions start\n')
parfor i = 1:length(variable_names)
        fprintf('i = %d\n',i) %Print current iteration number
    x = [param_input , results(:,i)]; %Construct matrix including the sampled parameter set 
                               %and the simulation results of each output
    [rho, pvalues] = partialcorr(x,'Type','Spearman'); %Calculate the PRCC
    rho_final(:,i) = rho(1:len_param, len_param+1); %Store the PRCC results 
    pvalues_final(:,i) = pvalues(1:len_param, len_param+1); %Store p-values
    fprintf('i = %d\n',i) %Print current iteration number
end

fprintf('PRCC solutions complete\n')
%Plot the PRCC heat map
figure(21)  
sensitivity_heatmap = heatmap(variable_names, parameter_names(metab_pars),rho_final(metab_pars,:),'colormap',jet);
 sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.CellLabelColor='none';
title('Global Sensitivity Analysis (PRCC): Metabolic Parameters')

figure(22)  
sensitivity_heatmap = heatmap(variable_names, parameter_names(cal_pars),rho_final(cal_pars,:),'colormap',jet);
 sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.CellLabelColor='none';
title('Global Sensitivity Analysis (PRCC): Calcium Parameters')
figure(23)  
sensitivity_heatmap = heatmap(variable_names, parameter_names(i_pars),rho_final(i_pars,:),'colormap',jet);
 sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.CellLabelColor='none';
title('Global Sensitivity Analysis (PRCC): INS Parameters')




%Plot the p-value heat map
figure(31)
sensitivity_heatmap = heatmap(variable_names, parameter_names(metab_pars),pvalues_final(metab_pars,:),'colormap',jet);
title('p-value Heat Map: Metabolic Parameters')
 sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.CellLabelColor='none';
figure(32)
sensitivity_heatmap = heatmap(variable_names, parameter_names(cal_pars),pvalues_final(cal_pars,:),'colormap',jet);
title('p-value Heat Map: Calcium Parameters')
 sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.CellLabelColor='none';
figure(33)
sensitivity_heatmap = heatmap(variable_names, parameter_names(i_pars),pvalues_final(i_pars,:),'colormap',jet);
title('p-value Heat Map: INS Parameters')
 sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.CellLabelColor='none';

 
%Plot PRCC heat map with only p < 0.05
rho_adjusted = rho_final; %copy PRCC matrix
for i = find(pvalues_final >= 0.05) %Set PRCC with p >= 0.05 to NaN
    rho_adjusted(i) = NaN;
end

figure(41)
sensitivity_heatmap = heatmap(variable_names, parameter_names(metab_pars),rho_adjusted(metab_pars,:),'colormap',jet);
title('PRCC Heat Map with p < 0.05: Metabolic Parameters')
sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.CellLabelColor='none';
figure(42)
sensitivity_heatmap = heatmap(variable_names, parameter_names(cal_pars),rho_adjusted(cal_pars,:),'colormap',jet);
title('PRCC Heat Map with p < 0.05: Calcium Parameters')
sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.CellLabelColor='none';
figure(43)
sensitivity_heatmap = heatmap(variable_names, parameter_names(i_pars),rho_adjusted(i_pars,:),'colormap',jet);
title('PRCC Heat Map with p < 0.05: INS Parameters')
sensitivity_heatmap.GridVisible = 'off';
sensitivity_heatmap.CellLabelColor='none';

%Plot 7D parallel plot 
figure(5)
idx = randi(N, 1, 150); %Randomly sample 150 simulation results out of N
parallelcoords(results(idx, [1,2,3,4]) , 'labels', {'pAMPK','pAKT','pmTORC1','pmTORC2'})
ylabel('Abundance')

toc %end timing
save globalSAresults %Store simulation results