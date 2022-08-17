%% Bifurcation Analysis
%Plot 3D bifurcation diagrams varying 2 parameters
%% Steps




%change IC of DEPTOR if needed
clc ; clear 
tic % start timer
getplotid
%Set initial parameters and range of bifurcation analysis
tt = 0:1000; %Timespanp
N = 1001; %Number of evalution points
Abundance = linspace(0.0001,1, N);%Re-define range of Abundance
%DEPTOR_IC = 1000; %Modify initial condition of DEPTOR
%y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;DEPTOR_IC;0;250;0;0;250;0]; %Initial condition
y0 = leung_icsource;
%Import model parameters
 load('param.mat');
pulsetime = 0;
pars = [leung_pars_source';param];
options = odeset('RelTol',1e-5, 'AbsTol',1e-8,  ...
    'InitialStep',1e-5, 'MaxOrder',5, 'BDF','on');
%Set points for V_IR in the range 1e-6 to 0.2 for simulation
param_1 = linspace(0.01, 1, 20);
param_2 = linspace(1, 10, 20);
param_3 = linspace(10, 30, 20);
param_range = [param_1(1:end-1),param_2(1:end-1), param_3];
kHYD = pars(4);
%Initialize matrices to store resulting steady state concentrations 
output_max_val = zeros(N, length(param_range), 2);
output_min_val = zeros(N, length(param_range), 2);
output_sum = zeros(N, length(param_range), 2);
output_product = zeros(N, length(param_range), 2);
output_div = zeros(N, length(param_range), 2);
output_pAMPK = zeros(N, length(param_range), 2);



%Initialize matrices to store resulting steady state concentrations 

%Solve the system ODEs at each V_IR and IC of [AMPK]
for j = 1:length(param_range)
    pars(105) = param_range(j); %C
    parfor i = 1:length(Abundance)
        y0_copy = y0; %Copy IC for parallel computing

        kHYD = Abundance(i); %Change kHYD
        
        [t, answer] = ode15s( @leung2022_replace,tt,y0_copy, options,pars, kHYD,pulsetime);%Solve the system of ODEs  

        %Find the max/min, sum, product, quotient of pmTORC1 and pmTORC2
        output_max_val(i, j, :) = max(answer(end-200:end, [pmTORC1, pmTORC2])); %Max of steady state
        output_min_val(i, j, :) = min(answer(end-200:end, [pmTORC1, pmTORC2])); %Min of steady state
		
%         %Max/min of steady state for pmTORC1 + pmTORC2
%         output_sum(i, j, :) = [max(answer(end-200:end, pmTORC1) + answer(end-200:end, pmTORC2)); 
%                                min(answer(end-200:end, pmTORC1) + answer(end-200:end, pmTORC2))];
% 		
%         %Max/min of steady state for pmTORC1 * pmTORC2                   
%         output_product(i, j, :) = [max(answer(end-200:end, pmTORC1) .* answer(end-200:end, pmTORC2));		
%                                   min(answer(end-200:end, pmTORC1) .* answer(end-200:end, pmTORC2))];		
% 		
%         %Max/min of steady state for pmTORC1 / pmTORC2                   
%         output_div(i, j, :) = [max(answer(end-200:end, pmTORC1) ./ answer(end-200:end, pmTORC2));		
%                               min(answer(end-200:end, pmTORC1) ./ answer(end-200:end, pmTORC2))];		
        
        %Max/min of steady state for pAMPK                  
        output_pAMPK(i, j, :) = max(answer(end-200:end, pAMPK   )); %Max of steady state
        output_pAMPK(i, j, :) = min(answer(end-200:end, pAMPK)); %Min of steady state
        
        fprintf('i = %d j = %d\n',i, j); % Print current indices of loop
    end
end


elapsedTime = toc % end timer
save(append('bifVIR_Khyd_', '1'))%Save simulation results