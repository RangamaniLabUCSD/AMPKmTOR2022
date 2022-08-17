%Local Sensitivity Analysis for the mTORC1/AMPK/ULK1 Pathway
clc; clear;            
tic %start timing

%Import model parameters


%Define parameters names for constructing the heatmap
parameter_names = {'KADP';'VmaxOP';'nH';'krest';'kstim';'kpost';'vAK';'kmt';'kmd';'kmm';'keqadk';'vmax20';...
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
'K_SIRT1';'K_SIRT1_by_pAMPK';'Km_SIRT1';'K_SIRT1_diss';'K_ULK1';'K_ULK1_by_pAMPK';'Km_ULK1';'K_pULK1';'K_pULK1_by_pmTORC1';'Km_pULK1';}; 
%parameter_select= [1:1:34,36:1:56,105:1:162];
%parameter_select = [1:1:3,7:1:10,31:1:34,36,37,38,40:1:55,105:162];
parameter_select = [1:1:3,31:1:34,36,37,38,105:162];
parameter_select = [1:162];

parameter_names = parameter_names(parameter_select);
parameter_names = strrep(parameter_names, '_', '\_'); % Replace '_' by '\_'

%Define variable names  
variable_names = {'Ca_C','pAMPK','pAKT','pmTORC1','pmTORC2'};  
%variable_names = {
%    'ATP','pAMPK','Ca_C','CAMKK2','pAKT','pmTORC1','pmTORC2','SIRT1','pULK1'};  
%Define the initial condition and time span for solving the system of ODEs
tt = 0:500; %Timespan
init = leung_ic ;
 load('param.mat');
pulsetime = 0;
options = odeset('RelTol',1e-5, 'AbsTol',1e-8,  ...
    'InitialStep',1e-5, 'MaxOrder',5, 'BDF','on');
pars = [leung_pars';param];

 kHYD = pars(4);
%Calculate the local sensitivty for each variable
Sensitivity = zeros(length(parameter_names), length(variable_names)); %Define matrix to store the sensitivities

%Solve the system of ODEs using the initial parameter set
[t, answer] = ode15s( @leung2022,tt,init, options,pars, kHYD,pulsetime);
S_p0_full = answer(end, :); % S(p0)

%Loop through each variable and parameter and find the sensitivty by
%approximating the derivatives and store the results in the sensitivity matrix

target_variables = [14,7,49,51,53];
target_params    = parameter_select;





for i = 1:length(target_variables)
    parfor j = 1:length(target_params)    
        
        %Set the vector param_modified with the modified parameter
        param_modified = pars; %copy original parameter set
        p = param_modified(target_params(j)); %Find the parameter of interest
        dp = 0.05 * p; %Set dp as 0.005*p
        param_modified(target_params(j)) = p + dp; %Perturb parameter of interest
        
        %Get S(p0)
        S_p0_copy = S_p0_full; %Copy for parallel computing 
        S_p0 = S_p0_copy(target_variables(i)); %Find S(p0)
       
        %Solve the system of ODEs with the modified parameters,and find the
        %terminal point of the solution for the variable of interest
[t, answer] = ode15s(@leung2022_replace,tt,init, options,param_modified, kHYD,pulsetime);
        S_p0_dp = answer(end, target_variables(i)); % Get S(p0+dp)
        
        % Calculate the relative local sensitivity
        % Sensitivity = p/S * dS/dp
        % where dS/dp approximated by dS/dp = [S(p0+dp)-S(p0)]/dp
        Sensitivity(j, i) = p/S_p0*(S_p0_dp - S_p0)/dp; 
        
        %Print current i and j
        fprintf('i = %d, j = %d\n', target_variables(i),target_params(j))
    end
end
 

%Construct the heat map
figure(1)
Sensitivity = abs(Sensitivity);
sensitivity_heatmap = heatmap(variable_names, parameter_select, Sensitivity, 'Colormap', parula);
sensitivity_heatmap.Title = 'Local Sensitivity Analysis';
% sensitivity_heatmap.GridVisible = 'off'; %Turn grid off
sensitivity_heatmap.FontSize = 12;
sensitivity_heatmap.ColorScaling ='log';

toc %stop timing

