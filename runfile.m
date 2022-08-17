
% main runfile AL
% 10/08/2021
% 

% clear workspace
clear
close all

load('param.mat')
%% Which Figures do you want to generate? 
% From 2-7
% For figure 8, need to use bifurcation scripts
% for figure 9, enter 9
figselect = [2,3,4,5,6,7,9,10];
getplotid


 
 %% Initial conditions
 init = leung_ic;
 pars = [leung_pars';param];

%Integration options
options = odeset('InitialStep',1e-6);
% options = odeset('RelTol',1e-5, 'AbsTol',1e-8, 'MaxStep',1e-1, ...
%     'InitialStep',1e-5, 'MaxOrder',5, 'BDF','on','NonNegative',[1:37,39:60]);
t_prior = 0 ;
t_start = 5000;
pulsetime = -100;


probe = 56;

init_baseval = init;


scan = 1.5;
    init(probe) = init_baseval(probe)*scan;
kHYD = pars(4);
        [t0, y_ss]= ode15s(@leung2022,[t_prior t_start],init,...
              options,pars, kHYD,pulsetime);

for fig_index = figselect
%% Initialization of Runtime Parameters
t_prior = 0 ;
t_start = 5000;
pulsetime = 0;
kHYD = pars(4);
glutconc = 100; % muM
disp(['Simulating for Figure ' num2str(fig_index)])


   %% Figure 2 : Steady State System

    if fig_index == 2  
    
        %% Simulation Engine for Steady State
        % initialization run
        % or load premade initialization data
%                 initcond =  y_ss(end,:);
%                 
%         %% Sensitivity Analysis
%         [t, y]= ode15s(@leung2022,[0,500],initcond,...
%               options,co_pars, kHYD,pulsetime);
t = t0;
y = y_ss;
        plotfigure2
    end
   %% Figure 3 Simple stimulus

    if fig_index == 3  
    
        %% Simulation Engine
        % initialization run
        % or load premade initialization data

        
        initcond =  y_ss(end,:);
        stimstart =0;


        y = [y_ss(end,:)];
        t = [0];
        stim_dur =50;
        pulsetime = 0;    
        frequency = 1; %Hz
        intervals = 1/frequency;%s 
        timevec   =  (stimstart:intervals:stimstart+stim_dur);
            
        for i = 1:(size(timevec,2)-1)        %simulate stimulus section     
            startpoint = timevec(i)+0.001;
            endpoint = timevec(i+1);
            initcond =  y(end,:);
            initcond(12) = glutconc;
            pulsetime = startpoint;
            [t_temp,x_temp]= ode15s(@leung2022,[startpoint endpoint],initcond,...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
        end
        recov_dur = 5000;
        recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
        [t_temp,x_temp]= ode15s(@leung2022,[recovtime],y(end,:),...
        options,pars, kHYD,pulsetime); % simulate steady state
        y = [y;x_temp];
        t = [t;t_temp];
        stimstart=recovtime(end)+0.001;
        t_auc = t;
        y_auc = y;
        
        t = [t0(end-500:end)-5000;t];
        y = [y_ss(end-500:end,:);y];
        %compute metrics
    
            [~,t200] = min(abs(t_auc-200));
            [~,t400] = min(abs(t_auc-400));
            [~,t600] = min(abs(t_auc-600));
            % compute AUC

            % compute baseline AUC
                AMPK_base = trapz(t_auc(t400:t600),y_auc(t400:t600,pAMPK));
                mTORC1_base = trapz(t_auc(t400:t600),y_auc(t400:t600,pmTORC1));
                mTORC2_base = trapz(t_auc(t400:t600),y_auc(t400:t600,pmTORC2));
                CaC_base = trapz(t_auc(t400:t600),y_auc(t400:t600,Ca_C));
            % compute AUC for time period
                AMPK_AUC = trapz(t_auc(1:t200),y_auc(1:t200,pAMPK));
                mTORC1_AUC = trapz(t_auc(1:t200),y_auc(1:t200,pmTORC1));
                mTORC2_AUC = trapz(t_auc(1:t200),y_auc(1:t200,pmTORC2));
                CaC_AUC = trapz(t_auc(1:t200),y_auc(1:t200,Ca_C));
%      % compute baseline AUC
%                 AMPK_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pAMPK));
%                 mTORC1_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pmTORC1));
%                 mTORC2_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pmTORC2));
%                 CaC_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),Ca_C));
%             % compute AUC for time period
%                 AMPK_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pAMPK));
%                 mTORC1_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pmTORC1));
%                 mTORC2_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pmTORC2));
%                 CaC_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),Ca_C));

 
        % normalize to baseline   
        normAMPK      = AMPK_AUC./AMPK_base    ;
        normmTOR1     = mTORC1_AUC./mTORC1_base    ;
        normmTOR2     = mTORC2_AUC./mTORC2_base    ;
        normCa        = CaC_AUC./CaC_base    ;
        t_TC = t_auc;        
        ampk_TC       = y_auc(:,pAMPK)  ;
        mTORC1_TC     = y_auc(:,pmTORC1)   ;
        mTORC2_TC     = y_auc(:,pmTORC2)   ;
        % time to peak

        

        
        %Period
        tmaxes_a = t_auc(islocalmax(ampk_TC,'MinSeparation',1));
        tmins_a  = t_auc(islocalmin(ampk_TC,'MinSeparation',1));
        tmaxes_m1 = t_auc(islocalmax(mTORC1_TC,'MinSeparation',1));
        tmins_m1  = t_auc(islocalmin(mTORC1_TC,'MinSeparation',1));     
        tmaxes_m2 = t_auc(islocalmax(mTORC2_TC,'MinSeparation',1));
        tmins_m2  = t_auc(islocalmin(mTORC2_TC,'MinSeparation',1));        
        
        dtmaxes_a =diff(tmaxes_a)    ;
        dtmaxes_m1=diff(tmaxes_m1)    ;
        dtmaxes_m2=diff(tmaxes_m2)    ;
        period_a= mean(dtmaxes_a(dtmaxes_a>3));   
        period_m1= mean(dtmaxes_m1(dtmaxes_m1>3));   
        period_m2 = mean(dtmaxes_m2(dtmaxes_m2>3));    
        %maximum amplitude
        %first compute mean values
        memaxval_a  =  mean(ampk_TC(islocalmax(ampk_TC,'MinSeparation',1)));  
        meminval_a  =  mean(ampk_TC(islocalmin(ampk_TC,'MinSeparation',1)));  
        memaxval_m1 =  mean(mTORC1_TC(islocalmax(mTORC1_TC,'MinSeparation',1))); 
        meminval_m1 =  mean(mTORC1_TC(islocalmin(mTORC1_TC,'MinSeparation',1))); 
        memaxval_m2 =  mean(mTORC2_TC(islocalmax(mTORC2_TC,'MinSeparation',1))); 
        meminval_m2 =  mean(mTORC2_TC(islocalmin(mTORC2_TC,'MinSeparation',1))); 

        maxamp_a= (max(ampk_TC(islocalmax(ampk_TC,'MinSeparation',1))) - mean([memaxval_a,meminval_a]))/mean([memaxval_a,meminval_a]);
        maxamp_m1= (max(mTORC1_TC(islocalmax(mTORC1_TC,'MinSeparation',1))) - mean([memaxval_m1,meminval_m1]))/mean([memaxval_m1,meminval_m1]);
        maxamp_m2= (max(mTORC2_TC(islocalmax(mTORC2_TC,'MinSeparation',1))) - mean([memaxval_m2,meminval_m2]))/mean([memaxval_m2,meminval_m2]);    

        %% Plots
        plotfigure3
        
    end      


   %% Figure 4 Frequency Scan



    if fig_index == 4  
    
    %% Simulation Engine
    % initialization run
    % or load premade initialization data
        %for freq_ind = 1:size(freq_vec,2)
        freq_scan = [0.1,1,10,100];
        
        for fscan_ind = 1:numel(freq_scan)
            j=freq_scan(fscan_ind);
            initcond =  y_ss(end,:);
            stimstart =0;
            y = [y_ss(end,:)];
            t = [0];
            stim_dur =5;        
            frequency = j; %Hz
            intervals = 1/frequency;%s 
            timevec   =  (stimstart:intervals:stimstart+stim_dur);
            if j==0.1
                timevec = [0,10];
            end
            
            for i = 1:(size(timevec,2)-1)             
                startpoint = timevec(i)+0.00001;
                endpoint = timevec(i+1);
                initcond =  y(end,:);
                initcond(12) = glutconc;
                pulsetime = startpoint;
                [t_temp,x_temp]= ode15s(@leung2022,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 5000;
            recovtime = timevec(end)+0.000001:0.1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;


        t = [t0(end-500:end)-5000;t];
        y = [y_ss(end-500:end,:);y];
        %% Plots
            
        plotfigure4

        end

    end    

   %% Figure 5

      if fig_index == 5  
    
    %% Simulation Engine
    % initialization run
    % or load premade initialization data
        %for freq_ind = 1:size(freq_vec,2)
        freq_scan = [0.1,1,10,100];
        
        for fscan_ind = 1:numel(freq_scan)
         j=freq_scan(fscan_ind);
            initcond =  y_ss(end,:);
            stimstart =0;
            y = [y_ss(end,:)];
            t = [0];
            stim_dur =5;        
            frequency = j; %Hz
            intervals = 1/frequency;%s 
            timevec   =  (stimstart:intervals:stimstart+stim_dur);
            if j==0.1
                timevec = [0,10];
            end
            
            for i = 1:(size(timevec,2)-1)             
                startpoint = timevec(i)+0.00001;
                endpoint = timevec(i+1);
                initcond =  y(end,:);
                initcond(12) = glutconc;
                pulsetime = startpoint;
                [t_temp,x_temp]= ode15s(@leung2022,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 5000;
            recovtime = timevec(end)+1e-9:0.1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;
            % find indexes for simulation time
            [~,t100] = min(abs(t-100));
            [~,t200] = min(abs(t-200));
            [~,t400] = min(abs(t-400));
            [~,t500] = min(abs(t-500));
            [~,t600] = min(abs(t-600));
            [~,t800] = min(abs(t-800));
            % compute AUC
            inttime = [t100,t200,t400];
            basetime = [t500,t600,t800];

            for inttime_ind = 1:numel(inttime)
            % compute baseline AUC
                AMPK_base(inttime_ind,fscan_ind) = trapz(t(t400:basetime(inttime_ind)),y(t400:basetime(inttime_ind),pAMPK));
                mTORC1_base(inttime_ind,fscan_ind) = trapz(t(t400:basetime(inttime_ind)),y(t400:basetime(inttime_ind),pmTORC1));
                mTORC2_base(inttime_ind,fscan_ind) = trapz(t(t400:basetime(inttime_ind)),y(t400:basetime(inttime_ind),pmTORC2));
                CaC_base(inttime_ind,fscan_ind) = trapz(t(t400:basetime(inttime_ind)),y(t400:basetime(inttime_ind),Ca_C));
            % compute AUC for time period
                AMPK_AUC(inttime_ind,fscan_ind) = trapz(t(1:inttime(inttime_ind)),y(1:inttime(inttime_ind),pAMPK));
                mTORC1_AUC(inttime_ind,fscan_ind) = trapz(t(1:inttime(inttime_ind)),y(1:inttime(inttime_ind),pmTORC1));
                mTORC2_AUC(inttime_ind,fscan_ind) = trapz(t(1:inttime(inttime_ind)),y(1:inttime(inttime_ind),pmTORC2));
                CaC_AUC(inttime_ind,fscan_ind) = trapz(t(1:inttime(inttime_ind)),y(1:inttime(inttime_ind),Ca_C));
%      % compute baseline AUC
%                 AMPK_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pAMPK));
%                 mTORC1_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pmTORC1));
%                 mTORC2_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pmTORC2));
%                 CaC_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),Ca_C));
%             % compute AUC for time period
%                 AMPK_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pAMPK));
%                 mTORC1_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pmTORC1));
%                 mTORC2_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pmTORC2));
%                 CaC_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),Ca_C));

            end    
        % normalize to baseline   
        normAMPK      = AMPK_AUC./AMPK_base    ;
        normmTOR1     = mTORC1_AUC./mTORC1_base    ;
        normmTOR2     = mTORC2_AUC./mTORC2_base    ;
        normCa        = CaC_AUC./CaC_base    ;
        t_TC = t;        
        ampk_TC       = y(:,pAMPK)  ;
        mTORC1_TC     = y(:,pmTORC1)   ;
        mTORC2_TC     = y(:,pmTORC2)   ;
        % time to peak

        

        
        %Period
        tmaxes_a = t(islocalmax(ampk_TC,'MinSeparation',100));
        tmins_a  = t(islocalmin(ampk_TC,'MinSeparation',100));
        tmaxes_m1 = t(islocalmax(mTORC1_TC,'MinSeparation',100));
        tmins_m1  = t(islocalmin(mTORC1_TC,'MinSeparation',100));     
        tmaxes_m2 = t(islocalmax(mTORC2_TC,'MinSeparation',100));
        tmins_m2  = t(islocalmin(mTORC2_TC,'MinSeparation',100));        
        
        dtmaxes_a =diff(tmaxes_a)    ;
        dtmaxes_m1=diff(tmaxes_m1)    ;
        dtmaxes_m2=diff(tmaxes_m2)    ;
        period_a(fscan_ind) = mean(dtmaxes_a(dtmaxes_a>3));   
        period_m1(fscan_ind) = mean(dtmaxes_m1(dtmaxes_m1>3));   
        period_m2(fscan_ind) = mean(dtmaxes_m2(dtmaxes_m2>3));    
        %maximum amplitude
        %first compute mean values
        memaxval_a  =  mean(ampk_TC(islocalmax(ampk_TC,'MinSeparation',100)));  
        meminval_a  =  mean(ampk_TC(islocalmin(ampk_TC,'MinSeparation',100)));  
        memaxval_m1 =  mean(mTORC1_TC(islocalmax(mTORC1_TC,'MinSeparation',100))); 
        meminval_m1 =  mean(mTORC1_TC(islocalmin(mTORC1_TC,'MinSeparation',100))); 
        memaxval_m2 =  mean(mTORC2_TC(islocalmax(mTORC2_TC,'MinSeparation',100))); 
        meminval_m2 =  mean(mTORC2_TC(islocalmin(mTORC2_TC,'MinSeparation',100))); 

        maxamp_a(fscan_ind)= (max(ampk_TC(islocalmax(ampk_TC,'MinSeparation',100))) - mean([memaxval_a,meminval_a]))/mean([memaxval_a,meminval_a]);
        maxamp_m1(fscan_ind)= (max(mTORC1_TC(islocalmax(mTORC1_TC,'MinSeparation',100))) - mean([memaxval_m1,meminval_m1]))/mean([memaxval_m1,meminval_m1]);
        maxamp_m2(fscan_ind)= (max(mTORC2_TC(islocalmax(mTORC2_TC,'MinSeparation',100))) - mean([memaxval_m2,meminval_m2]))/mean([memaxval_m2,meminval_m2]);

        %
        
        
        
        
        %% Plots
            
%        plotfigure4
        %figure(411)
        %plotfigure4_fft

        end
        % barlabel = categorical({'Calcium','pAMPK','pmTORC1','pmTORC2'});
        %bar(barlabel,[CaC_AUC./max(CaC_AUC);AMPK_AUC./max(AMPK_AUC);mTORC1_AUC./max(mTORC1_AUC);mTORC2_AUC./max(mTORC2_AUC)])
        %ylabel('Normalized Area Under Curve')
        %legend('0.1 Hz','1','10','100')
        plotfigure5_metrics
    end    

   %% Figure 6

   if fig_index == 6  
    clear
        %% Simulation Engine
        % initialization run
        % or load premade initialization data
        %for freq_ind = 1:size(freq_vec,2)
        parscan = [0.1,1,2,3];
        kHYDbasal = kHYD;
        for j = 1:numel(parscan)
            h = parscan(j);
            kHYD = kHYDbasal * h;
            initcond =  y_ss(end,:);
            stimstart =0;
            y = [y_ss(end,:)];
            t = [0];
            stim_dur =5;        
                frequency = 10; %Hz
                intervals = 1/frequency;%s 
                timevec   =  (stimstart:intervals:stimstart+stim_dur);
                
                for i = 1:(size(timevec,2)-1)             
                    startpoint = timevec(i)+0.001;
                    endpoint = timevec(i+1);
                    initcond =  y(end,:);
                    initcond(12) = glutconc;
                    pulsetime = startpoint;
                    [t_temp,x_temp]= ode15s(@leung2022,[startpoint endpoint],initcond,...
                    options,pars, kHYD,pulsetime);
                    y = [y;x_temp];
                    t = [t;t_temp];
                end
                recov_dur = 300;
                recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
                [t_temp,x_temp]= ode15s(@leung2022,[recovtime],y(end,:),...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
                stimstart=recovtime(end)+0.001;
            
            recov_dur = 300;
            recovtime = recovtime(end)+1:recovtime(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];        
            t_auc = t;
            y_auc = y;
            
            t = [t0(end-500:end)-5000;t];
        y = [y_ss(end-500:end,:);y];
        %% Plots


            [~,t200] = min(abs(t_auc-200));
            [~,t400] = min(abs(t_auc-400));
            [~,t600] = min(abs(t_auc-600));
            % compute AUC

            % compute baseline AUC
                AMPK_base(j) = trapz(t_auc(t400:t600),y_auc(t400:t600,pAMPK));
                mTORC1_base(j) = trapz(t_auc(t400:t600),y_auc(t400:t600,pmTORC1));
                mTORC2_base(j) = trapz(t_auc(t400:t600),y_auc(t400:t600,pmTORC2));
                CaC_base(j) = trapz(t_auc(t400:t600),y_auc(t400:t600,Ca_C));
            % compute AUC for time period
                AMPK_AUC(j) = trapz(t_auc(1:t200),y_auc(1:t200,pAMPK))
                mTORC1_AUC(j) = trapz(t_auc(1:t200),y_auc(1:t200,pmTORC1));
                mTORC2_AUC(j) = trapz(t_auc(1:t200),y_auc(1:t200,pmTORC2));
                CaC_AUC(j) = trapz(t_auc(1:t200),y_auc(1:t200,Ca_C));
%      % compute baseline AUC
%                 AMPK_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pAMPK));
%                 mTORC1_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pmTORC1));
%                 mTORC2_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),pmTORC2));
%                 CaC_base(inttime_ind,fscan_ind) = trapz(t(t400:t600),y(t400:basetime(inttime_ind),Ca_C));
%             % compute AUC for time period
%                 AMPK_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pAMPK));
%                 mTORC1_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pmTORC1));
%                 mTORC2_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),pmTORC2));
%                 CaC_AUC(inttime_ind,fscan_ind) = trapz(t(1:t200),y(1:inttime(inttime_ind),Ca_C));

 
        % normalize to baseline   

        t_TC = t_auc;        
        ampk_TC       = y_auc(:,pAMPK)  ;
        mTORC1_TC     = y_auc(:,pmTORC1)   ;
        mTORC2_TC     = y_auc(:,pmTORC2)   ;
        % time to peak

        %steady states
        AMPK_SS(j)   =  mean(ampk_TC(end-10:end));
        mTOR1_SS(j)  =  mean(mTORC1_TC(end-10:end));
        mTOR2_SS(j)  =  mean(mTORC2_TC(end-10:end));

        
        %Period
         %subset data
         tframe = t200:5753;
         t_stab = t_auc(tframe);
        AMPK_stab  =  ampk_TC(tframe)./mean(ampk_TC(tframe));
        mTOR1_stab  =  mTORC1_TC(tframe)./mean(mTORC1_TC(tframe));
        mTOR2_stab  =  mTORC2_TC(tframe)./mean(mTORC2_TC(tframe));
            % remove by average
            [apeaks,alocs] = findpeaks(AMPK_stab);
            [m1peaks,m1locs] = findpeaks(mTOR1_stab);
            [m2peaks,m2locs] = findpeaks(mTOR2_stab);


            [apeaks,alocs] = findpeaks(AMPK_stab,'MinPeakHeight',1+0.1*(max(apeaks)-1));
            [m1peaks,m1locs] = findpeaks(mTOR1_stab,'MinPeakHeight',1+0.1*(max(m1peaks)-1));
            [m2peaks,m2locs] = findpeaks(mTOR2_stab,'MinPeakHeight',1+0.1*(max(m2peaks)-1));
            %determine peak height
            
            %peaks by minimum distance

            %return mean period
        period_a(j)=  mean(diff(t_auc(tframe(alocs))));   
        period_m1(j)=  mean(diff(t_auc(tframe(m1locs))));   
        period_m2(j) = mean(diff(t_auc(tframe(m2locs))));    
        %maximum amplitude
        %first compute mean values
%         memaxval_a  =  mean(ampk_TC(islocalmax(ampk_TC,'MinSeparation',1)));  
%         meminval_a  =  mean(ampk_TC(islocalmin(ampk_TC,'MinSeparation',1)));  
%         memaxval_m1 =  mean(mTORC1_TC(islocalmax(mTORC1_TC,'MinSeparation',1))); 
%         meminval_m1 =  mean(mTORC1_TC(islocalmin(mTORC1_TC,'MinSeparation',1))); 
%         memaxval_m2 =  mean(mTORC2_TC(islocalmax(mTORC2_TC,'MinSeparation',1))); 
%         meminval_m2 =  mean(mTORC2_TC(islocalmin(mTORC2_TC,'MinSeparation',1))); 

        maxamp_a(j)= (max(ampk_TC(islocalmax(ampk_TC,'MinSeparation',1))) );
        maxamp_m1(j)= (min(mTORC1_TC(islocalmin(mTORC1_TC,'MinSeparation',1))) );
        maxamp_m2(j)= (max(mTORC2_TC(islocalmax(mTORC2_TC,'MinSeparation',1))) );    

        plotfigure6

        end
                AMPK_SS  = AMPK_SS./AMPK_SS(2)   ;
                mTOR1_SS  = mTOR1_SS./mTOR1_SS(2)   ;
                mTOR2_SS  = mTOR2_SS./mTOR2_SS(2)   ;
               maxamp_a   = maxamp_a./maxamp_a(2)/100;
                maxamp_m1 = maxamp_m1./maxamp_m1(2)/100;
               maxamp_m2  =    maxamp_m2./maxamp_m2(2)/100;
                

        normAMPK      = AMPK_AUC./AMPK_AUC(2)    ;
        normmTOR1     = mTORC1_AUC./mTORC1_AUC(2) ;   
        normmTOR2     = mTORC2_AUC./mTORC2_AUC(2)  ;  
        normCa        = CaC_AUC./CaC_base(2)    ;


                plotfigure6_metrics

    end    

%% Figure 7: VIR
   if fig_index == 7  
        clear AMPK_AUC mTORC1_AUC CaC_AUC mTORC2_AUC
        %% Simulation Engine
        % initialization run
        % or load premade initialization data
        [~,y_ss]= ode15s(@leung2022_ampkact,[0:1:1000],init,options,pars, kHYD,pulsetime,0.1,0.1);
            initcond =  y_ss(end,:);
         [tbaseline,ybaseline]= ode15s(@leung2022_ampkact,[0:1:200],init,options,pars, kHYD,pulsetime,0.1,0.1);
            AMPK_base =  trapz(tbaseline,ybaseline);
        mTORC1_base=  trapz(tbaseline,ybaseline);
        mTORC2_base =  trapz(tbaseline,ybaseline);
        CaC_base=  trapz(tbaseline,ybaseline);
        %for freq_ind = 1:size(freq_vec,2)
        parscan = [0.1,1,2,10];
        khyd_vary = [0.1,1,2,10];
        colorvector = ['#0072BD';'#D95319';'#EDB120';	'#7E2F8E'];
        markervector = {'-';'--';':';'-.'};      
        kHYDbasal = kHYD;
for VIR_ind = 1:numel(parscan)
        for khyd_ind = 1:numel(khyd_vary)    
            h = khyd_vary(khyd_ind);
            kHYD = kHYDbasal * h;
            freepar = parscan(VIR_ind);
            freepar_hyd = khyd_vary(khyd_ind);
            [~,y_ss]= ode15s(@leung2022_ampkact,[0:1:1000],init,options,pars, kHYD,pulsetime,freepar,freepar_hyd);
            initcond =  y_ss(end,:);
            stimstart =0;
            y = [y_ss(end,:)];
            t = [0];
            stim_dur =5;

                
                frequency = 10; %Hz
                intervals = 1/frequency;%s 
                timevec   =  (stimstart:intervals:stimstart+stim_dur);
                
                for i = 1:(size(timevec,2)-1)             
                    startpoint = timevec(i)+0.001;
                    endpoint = timevec(i+1);
                    initcond =  y(end,:);
                    initcond(12) = glutconc;
                    pulsetime = startpoint;
                    [t_temp,x_temp]= ode15s(@leung2022_ampkact,[startpoint endpoint],initcond,...
                    options,pars, kHYD,pulsetime,freepar,freepar_hyd);
                    y = [y;x_temp];
                    t = [t;t_temp];
                end
                recov_dur = 500;
                recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
                [t_temp,x_temp]= ode15s(@leung2022_ampkact,[recovtime],y(end,:),...
                options,pars, kHYD,pulsetime,freepar,freepar_hyd);
                y = [y;x_temp];
                t = [t;t_temp];
                stimstart=recovtime(end)+0.001;
            
            recov_dur = 500;
            recovtime = recovtime(end)+1:recovtime(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_ampkact,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime,freepar,freepar_hyd);
            y = [y;x_temp];
            t = [t;t_temp];
            figure(8)
            hold on
        [~,t200] = min(abs(t-200));
        [~,t400] = min(abs(t-400));
        [~,t600] = min(abs(t-600));
       

        
        AMPK_AUC(VIR_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pAMPK));
        mTORC1_AUC(VIR_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pmTORC1));
        mTORC2_AUC(VIR_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pmTORC2));
        CaC_AUC(VIR_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,Ca_C));

%         AMPK_AUC = AMPK_AUC ./ AMPK_base;
%         mTORC1_AUC = mTORC1_AUC ./ mTORC1_base;
%         mTORC2_AUC = mTORC2_AUC ./ mTORC2_base;
%         CaC_AUC = CaC_AUC ./ CaC_base;
        
            t = [t0(end-500:end)-5000;t];
        y = [y_ss(end-500:end,:);y];
                plotfigure7
        %% Plots

        end
        end

                
    end    

if (fig_index == 9 )
    % Frequency vs VIR
    clear AMPK_AUC mTORC1_AUC CaC_AUC mTORC2_AUC
    %% Simulation Engine
    % initialization run
    % or load premade initialization data
    [~,y_ss]= ode15s(@leung2022_ampkact,[0:1:1000],init,options,pars, kHYD,pulsetime,0.1,0.1);
    initcond =  y_ss(end,:);
    [tbaseline,ybaseline]= ode15s(@leung2022_ampkact,[0:1:200],init,options,pars, kHYD,pulsetime,0.1,0.1);
    AMPK_base =  trapz(tbaseline,ybaseline);
    mTORC1_base=  trapz(tbaseline,ybaseline);
    mTORC2_base =  trapz(tbaseline,ybaseline);
    CaC_base=  trapz(tbaseline,ybaseline);
    %for freq_ind = 1:size(freq_vec,2)
    parscan = [0.1,1,2,10];
    freq_vary = [0.1,1,10,50];
    colorvector = ['#0072BD';'#D95319';'#EDB120';	'#7E2F8E'];
    markervector = {'-';'--';':';'-.'};      
    for VIR_ind = 1:numel(parscan)
        for freq_ind = 1:numel(freq_vary)    
            h = freq_vary(freq_ind);
            freepar = parscan(VIR_ind);
            [~,y_ss]= ode15s(@leung2022_ampkact,[0:1:1000],init,options,pars, kHYD,pulsetime,freepar,1);
            initcond =  y_ss(end,:);
            stimstart =0;
            y = [y_ss(end,:)];
            t = [0];
            stim_dur =5;
            
            
            frequency = h; %Hz
            intervals = 1/frequency;%s 
            timevec   =  (stimstart:intervals:stimstart+stim_dur);
            
            for i = 1:(size(timevec,2)-1)             
                startpoint = timevec(i)+0.001;
                endpoint = timevec(i+1);
                initcond =  y(end,:);
                initcond(12) = glutconc;
                pulsetime = startpoint;
                [t_temp,x_temp]= ode15s(@leung2022_ampkact,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime,freepar,1);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 500;
            recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_ampkact,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime,freepar,1);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;
            recov_dur = 500;
            recovtime = recovtime(end)+1:recovtime(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_ampkact,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime,freepar,1);
            y = [y;x_temp];
            t = [t;t_temp];
            [~,t200] = min(abs(t-200));
            [~,t400] = min(abs(t-400));
            [~,t600] = min(abs(t-600));
            
            
            
            AMPK_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,pAMPK));
            mTORC1_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,pmTORC1));
            mTORC2_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,pmTORC2));
            CaC_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,Ca_C));
            
            %         AMPK_AUC = AMPK_AUC ./ AMPK_base;
            %         mTORC1_AUC = mTORC1_AUC ./ mTORC1_base;
            %         mTORC2_AUC = mTORC2_AUC ./ mTORC2_base;
            %         CaC_AUC = CaC_AUC ./ CaC_base;
            
            t = [t0(end-500:end)-5000;t];
            y = [y_ss(end-500:end,:);y];
            plotfigure9
            %% Plots
        end
    end
end
if (fig_index == 10 )
    % Frequency vs kHYD
    clear AMPK_AUC mTORC1_AUC CaC_AUC mTORC2_AUC
    %% Simulation Engine
    % initialization run
    % or load premade initialization data
    [~,y_ss]= ode15s(@leung2022_ampkact,[0:1:1000],init,options,pars, kHYD,pulsetime,0.1,0.1);
    initcond =  y_ss(end,:);
    [tbaseline,ybaseline]= ode15s(@leung2022_ampkact,[0:1:200],init,options,pars, kHYD,pulsetime,0.1,0.1);
    AMPK_base =  trapz(tbaseline,ybaseline);
    mTORC1_base=  trapz(tbaseline,ybaseline);
    mTORC2_base =  trapz(tbaseline,ybaseline);
    CaC_base=  trapz(tbaseline,ybaseline);
    %for freq_ind = 1:size(freq_vec,2)
    freq_vary = [0.1,1,10,50];
    khyd_vary = [0.1,1,2,10];
    colorvector = ['#0072BD';'#D95319';'#EDB120';	'#7E2F8E'];
    markervector = {'-';'--';':';'-.'};      
    kHYDbasal = kHYD;
    for freq_ind = 1:numel(freq_vary)
        for khyd_ind = 1:numel(khyd_vary)    
            
            h = khyd_vary(freq_ind);
            kHYD = kHYDbasal * h;
            frequency = freq_vary(freq_ind);
            freepar_hyd = khyd_vary(khyd_ind);
            [~,y_ss]= ode15s(@leung2022_ampkact,[0:1:1000],init,options,pars, kHYD,pulsetime,1,freepar_hyd);
            initcond =  y_ss(end,:);
            stimstart =0;
            y = [y_ss(end,:)];
            t = [0];
            stim_dur =5;               
            intervals = 1/frequency;%s 
            timevec   =  (stimstart:intervals:stimstart+stim_dur);
            
            for i = 1:(size(timevec,2)-1)             
                startpoint = timevec(i)+0.001;
                endpoint = timevec(i+1);
                initcond =  y(end,:);
                initcond(12) = glutconc;
                pulsetime = startpoint;
                [t_temp,x_temp]= ode15s(@leung2022_ampkact,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime,1,freepar_hyd);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 500;
            recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_ampkact,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime,1,freepar_hyd);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;
            
            recov_dur = 500;
            recovtime = recovtime(end)+1:recovtime(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_ampkact,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime,1,freepar_hyd);
            y = [y;x_temp];
            t = [t;t_temp];

            hold on
            [~,t200] = min(abs(t-50));
            [~,t400] = min(abs(t-400));
            [~,t600] = min(abs(t-600));
            
            
            
            AMPK_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pAMPK));
            mTORC1_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pmTORC1));
            mTORC2_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pmTORC2));
            CaC_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,Ca_C));
            
            %         AMPK_AUC = AMPK_AUC ./ AMPK_base;
            %         mTORC1_AUC = mTORC1_AUC ./ mTORC1_base;
            %         mTORC2_AUC = mTORC2_AUC ./ mTORC2_base;
            %         CaC_AUC = CaC_AUC ./ CaC_base;
            
            t = [t0(end-500:end)-5000;t];
            y = [y_ss(end-500:end,:);y];
            plotfigure10
            %% Plots
        end
    end
end
 %       legend('1','2','3','4','5','6','7')
end