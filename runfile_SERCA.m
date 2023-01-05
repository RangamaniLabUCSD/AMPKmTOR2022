
% main runfile AL
% 10/08/2021
% 

% clear workspace
clear
% close all

load('finalpars.mat')
load('newparams.mat')
%% Which Figures do you want to generate? 
% From 2-7
% For figure 8, need to use bifurcation scripts

figselect = [9];
figselect = [2,3,4,5,6,7,9,10];


% figselect = [6];
getplotid
% %figselect = [2];
 finalpars(105) = 3;
pars = finalpars;

 %% Initial conditions
 init = leung_icsource;
 
 %% inserted to swap to fit pars
%  pars = newparams(1:162);
%  finalpars =pars;
% init = newparams(163:222);
% init(38) =-init(38);
%Integration options
options = odeset('InitialStep',1e-6);
% options = odeset('RelTol',1e-5, 'AbsTol',1e-8, 'MaxStep',1e-1, ...
%     'InitialStep',1e-5, 'MaxOrder',5, 'BDF','on','NonNegative',[1:37,39:60]);
t_prior = 0 ;
t_start = 10000;
pulsetime = -100;
init_baseval = init;




kHYD = pars(4);
        [t0, y_ss]= ode15s(@leung2022_SERCA,[t_prior t_start],init,...
              options,pars, kHYD,pulsetime);

for fig_index = figselect
%% Initialization of Runtime Parameters
pulsetime = 0;
kHYD = pars(4);
glutconc = 50; % muM
disp(['Simulating for Figure ' num2str(fig_index)])


   %% Figure 2 : Steady State System

    if fig_index == 2  
    
        %% Simulation Engine for Steady State
        % initialization run
        % or load premade initialization data
%                 initcond =  y_ss(end,:);
%                 
%         %% Sensitivity Analysis
%         [t, y]= ode15s(@leung2022_SERCA,[0,500],initcond,...
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
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
        end
        recov_dur = 5000;
        recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
        [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
        options,pars, kHYD,pulsetime); % simulate steady state
        y = [y;x_temp];
        t = [t;t_temp];
        stimstart=recovtime(end)+0.001;
        t_auc = t;
        y_auc = y;
        
        t = [t0(1:end)-t_start;t];
        y = [y_ss(1:end,:);y];
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
        ampk_based = ampk_TC - mean(ampk_TC(floor(7*end/8):end));
        m1_based = mTORC1_TC - mean(mTORC1_TC(floor(7*end/8):end));
        m2_based = mTORC2_TC - mean(mTORC2_TC(floor(7*end/8):end));
        [peakheight_a,tmaxes_a] = findpeaks(ampk_based,t_TC,'MinPeakWidth',2);
        [mins_a,tmins_a]  = findpeaks(-ampk_based,t_TC,'MinPeakWidth',2);
        [peakheight_m1,tmaxes_m1] =findpeaks(m1_based,t_TC,'MinPeakWidth',2);
        [ mins_m1,tmins_m1]  = findpeaks(-m1_based,t_TC,'MinPeakWidth',2);     
        [peakheight_m2,tmaxes_m2] = findpeaks(m2_based,t_TC,'MinPeakWidth',2);
        [mins_m2,tmins_m2]  = findpeaks(-m2_based,t_TC,'MinPeakWidth',2);        


        


        dtmaxes_a =diff(tmaxes_a)    ;
        dtmaxes_m1=diff(tmaxes_m1)    ;
        dtmaxes_m2=diff(tmaxes_m2)    ;
        period_a= mean(dtmaxes_a(dtmaxes_a>3));   
        period_m1= mean(dtmaxes_m1(dtmaxes_m1>3));   
        period_m2 = mean(dtmaxes_m2(dtmaxes_m2>3));    
        %maximum amplitude
        %first compute mean peak values

        maxamp_a=  max(peakheight_a);
        maxamp_m1= max(peakheight_m1) ;
        maxamp_m2= max(peakheight_m2) ;    
 
        %% Plots
        plotfigure3
        
    end      


   %% Figure 4 Frequency Scan



    if fig_index == 4  
    
    %% Simulation Engine
    % initialization run
    % or load premade initialization data
        %for freq_ind = 1:size(freq_vec,2)
      %  freq_scan = [0.1,1,10,100];
              freq_scan = [0.1,1,10,50];

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
                [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 5000;
            recovtime = timevec(end)+0.000001:0.1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;


        t = [t0(1:end)-t_start;t];
        y = [y_ss(1:end,:);y];
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
        freq_scan = [0.1,1,10,50];
        
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
                [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 5000;
            recovtime = timevec(end)+1e-9:0.1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;
            % find indexes for simulation time

            [~,t5] = min(abs(t-5));
            [~,t100] = min(abs(t-50));
            [~,t200] = min(abs(t-200));
            [~,t400] = min(abs(t-400));
            [~,t500] = min(abs(t-500));
            [~,t700] = min(abs(t-700));
            [~,t800] = min(abs(t-750));
            % compute AUC
            inttime = [t100,t200,t400];
            basetime = t800;


            % compute baseline AUC
                AMPK_base(1,fscan_ind) = trapz(t(t700:t800),y(t700:t800,pAMPK));
                mTORC1_base(1,fscan_ind) = trapz(t(t700:t800),y(t700:t800,pmTORC1));
                mTORC2_base(1,fscan_ind) = trapz(t(t700:t800),y(t700:t800,pmTORC2));
                CaC_base(1,fscan_ind) = trapz(t(t700:t800),y(t700:t800,Ca_C));
            % compute AUC for time period
                AMPK_AUC(1,fscan_ind) = trapz(t(1:t100),y(1:t100,pAMPK));
                mTORC1_AUC(1,fscan_ind) = trapz(t(1:t100),y(1:t100,pmTORC1));
                mTORC2_AUC(1,fscan_ind) = trapz(t(1:t100),y(1:t100,pmTORC2));
                CaC_AUC(1,fscan_ind) = trapz(t(1:t100),y(1:t100,Ca_C));

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

        

        %Decay time
        %Period
ampk_ss = mean(ampk_TC(floor(95*end/100):end));
m1_ss =mean(mTORC1_TC(floor(95*end/100):end));
m2_ss =mean(mTORC2_TC(floor(95*end/100):end));

ampk_based = ampk_TC - ampk_ss;
        m1_based = mTORC1_TC - m1_ss;
        m2_based = mTORC2_TC - m2_ss;
        [peakheight_a,tmaxes_a] = findpeaks(ampk_based,t_TC,'MinPeakWidth',2);
        [mins_a,tmins_a]  = findpeaks(-ampk_based,t_TC,'MinPeakWidth',2);
        [peakheight_m1,tmaxes_m1] =findpeaks(m1_based,t_TC,'MinPeakWidth',2);
        [ mins_m1,tmins_m1]  = findpeaks(-m1_based,t_TC,'MinPeakWidth',2);     
        [peakheight_m2,tmaxes_m2] = findpeaks(m2_based,t_TC,'MinPeakWidth',2);
        [mins_m2,tmins_m2]  = findpeaks(-m2_based,t_TC,'MinPeakWidth',2);  

        %Algorithm to find point when it timecourse returns to baseline
        %values

        %  take abs value of ss subtracted values
        ab_ampktc = abs(ampk_based);
        ab_m1tc = abs(m1_based);
        ab_m2tc = abs(m2_based);
        %  moving window mean of time course with a step size of 30 seconds
        windowstart = t5;

         ampk_movmean = movmean(ab_ampktc(windowstart:end), 60);
         m1_movmean = movmean(ab_m1tc(windowstart:end), 60);
         m2_movmean = movmean(ab_m2tc(windowstart:end), 60);
        % < some threshold value
        eq_thresh = 1e-8;
        % find closest peak value
roc_AMPK = abs(diff(ampk_TC)./diff(t_TC));
roc_m1 = abs(diff(mTORC1_TC)./diff(t_TC));
roc_m2 = abs(diff(mTORC2_TC)./diff(t_TC));
ampkdiff_movmean = movmean(roc_AMPK(windowstart:end), 10);
m1diff_movmean = movmean(roc_m1(windowstart:end), 10);
m2diff_movmean = movmean(roc_m2(windowstart:end), 10);

   threshind_a = find(ampkdiff_movmean-eq_thresh<0);
         threshind_m1 = find(m1diff_movmean-eq_thresh<0);
         threshind_m2 = find(m2diff_movmean-eq_thresh<0);

%         
         tte_a(fscan_ind)= t_TC (windowstart + threshind_a(1));
         tte_m1(fscan_ind)= t_TC (windowstart + threshind_m1(1)) ;
         tte_m2(fscan_ind) =t_TC (windowstart + threshind_m2(1))   ;

        dtmaxes_a =diff(tmaxes_a)    ;
        dtmaxes_m1=diff(tmaxes_m1)    ;
        dtmaxes_m2=diff(tmaxes_m2)    ;
  
        %maximum amplitude
        %first compute mean values


        maxamp_a(fscan_ind)= max(peakheight_a)/ampk_ss;
        maxamp_m1(fscan_ind)= max(peakheight_m1)/m1_ss;
        maxamp_m2(fscan_ind)= max(peakheight_m2)/m2_ss;

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
    
        clear AMPK_AUC mTORC1_AUC CaC_AUC mTORC2_AUC

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
                    [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
                    options,pars, kHYD,pulsetime);
                    y = [y;x_temp];
                    t = [t;t_temp];
                end
                recov_dur = 1500;
                recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
                [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];


            t_auc = t;
            y_auc = y;
            
             t = [t0(1:end)-t_start;t];
        y = [y_ss(1:end,:);y];
        %% Plots

            [~,t50] = min(abs(t_auc-50));

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
                AMPK_AUC(j) = trapz(t_auc(1:t200),y_auc(1:t200,pAMPK));
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
%              normAMPK      = AMPK_AUC./AMPK_base    ;
%         normmTOR1     = mTORC1_AUC./mTORC1_base    ;
%         normmTOR2     = mTORC2_AUC./mTORC2_base    ;
%         normCa        = CaC_AUC./CaC_base    ;
         
        ampk_TC       = y_auc(:,pAMPK)  ;
        mTORC1_TC     = y_auc(:,pmTORC1)   ;
        mTORC2_TC     = y_auc(:,pmTORC2)   ;
        % time to peak

        

        %Decay time
        %Period
ampk_ss = mean(ampk_TC(floor(90*end/100):end));
m1_ss =mean(mTORC1_TC(floor(90*end/100):end));
m2_ss =mean(mTORC2_TC(floor(90*end/100):end));

ampk_based = ampk_TC - ampk_ss;
        m1_based = mTORC1_TC - m1_ss;
        m2_based = mTORC2_TC - m2_ss;
        [peakheight_a,tmaxes_a] = findpeaks(ampk_based,t_TC,'MinPeakWidth',2);
        [mins_a,tmins_a]  = findpeaks(-ampk_based,t_TC,'MinPeakWidth',2);
        [peakheight_m1,tmaxes_m1] =findpeaks(m1_based,t_TC,'MinPeakWidth',2);
        [ mins_m1,tmins_m1]  = findpeaks(-m1_based,t_TC,'MinPeakWidth',2);     
        [peakheight_m2,tmaxes_m2] = findpeaks(m2_based,t_TC,'MinPeakWidth',2);
        [mins_m2,tmins_m2]  = findpeaks(-m2_based,t_TC,'MinPeakWidth',2);  
windowstart = t50;
roc_AMPK = abs(diff(ampk_TC)./diff(t_TC));
roc_m1 = abs(diff(mTORC1_TC)./diff(t_TC));
roc_m2 = abs(diff(mTORC2_TC)./diff(t_TC));
ampkdiff_movmean = movmean(roc_AMPK(windowstart:end), 10);
m1diff_movmean = movmean(roc_m1(windowstart:end), 10);
m2diff_movmean = movmean(roc_m2(windowstart:end), 10);
  eq_thresh = 1e-6;
   threshind_a = find(ampkdiff_movmean-eq_thresh<0);
         threshind_m1 = find(m1diff_movmean-eq_thresh<0);
         threshind_m2 = find(m2diff_movmean-eq_thresh<0);
%         %Algorithm to find point when it timecourse returns to baseline
%         %values
% 
%         %  take abs value of ss subtracted values
%         ab_ampktc = abs(ampk_based);
%         ab_m1tc = abs(m1_based);
%         ab_m2tc = abs(m2_based);
%         %  moving window mean of time course with a step size of 30 seconds
%         windowstart = t50;
%          ampk_movmean = movmean(ab_ampktc(windowstart:end), 30);
%          m1_movmean = movmean(ab_m1tc(windowstart:end), 30);
%          m2_movmean = movmean(ab_m2tc(windowstart:end), 30);
%         % < some threshold value
%         eq_thresh = 1e-6;
%         % find closest peak value
%         threshind_a = find(ampk_movmean-eq_thresh<0);
%         threshind_m1 = find(m1_movmean-eq_thresh<0);
%         threshind_m2 = find(m2_movmean-eq_thresh<0);

%         % put into period for quick plotting, change in plotfigure5_metrics
%         % later
%         
         tte_a(j)= t_TC (windowstart + threshind_a(1));
         tte_m1(j)= t_TC (windowstart + threshind_m1(1)) ;
         tte_m2(j) =t_TC (windowstart + threshind_m2(1))   ;
%         if isempty(tte_a)
%             tte_a(fscan_ind) = 1;
%         end




        dtmaxes_a =diff(tmaxes_a)    ;
        dtmaxes_m1=diff(tmaxes_m1)    ;
        dtmaxes_m2=diff(tmaxes_m2)    ;
%         period_a= mean(dtmaxes_a(dtmaxes_a>3));   
%         period_m1= mean(dtmaxes_m1(dtmaxes_m1>3));   
%         period_m2 = mean(dtmaxes_m2(dtmaxes_m2>3));    
        %maximum amplitude
        %first compute mean peak values


 

        
%             [freq_a,zeta_a] = damp(t_TC,ampk_TC);
%             tau_a = decayTime(zeta_a,freq_a);


   
%         period_a(fscan_ind) = mean(dtmaxes_a(dtmaxes_a>3));   
%         period_m1(fscan_ind) = mean(dtmaxes_m1(dtmaxes_m1>3));   
%         period_m2(fscan_ind) = mean(dtmaxes_m2(dtmaxes_m2>3));    
        %maximum amplitude
        %first compute mean values


        maxamp_a(j)= max(peakheight_a)/ampk_ss;
        maxamp_m1(j)= max(peakheight_m1)/m1_ss;
        maxamp_m2(j)= max(peakheight_m2)/m2_ss;

        plotfigure6
                AMPK_SS(j) = ampk_ss;
                mTOR1_SS(j)  = m1_ss;
                mTOR2_SS(j)  =m2_ss;
        end %end parscan
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
                parscan = [0.1,1,2,10];
        khyd_vary = [0.1,1,2,10];
pulsetime =0;
        %% Simulation Engine
        % initialization run
        % or load premade initialization data
        [~,y_ss]= ode15s(@leung2022_SERCA,[0:1:1000],init,options,pars, kHYD,pulsetime);
            initcond =  y_ss(end,:);
         [tbaseline,ybaseline]= ode15s(@leung2022_SERCA,[0:1:200],initcond,options,pars, kHYD,pulsetime);

            AMPK_base =  trapz(tbaseline,ybaseline);
        mTORC1_base=  trapz(tbaseline,ybaseline);
        mTORC2_base =  trapz(tbaseline,ybaseline);
        CaC_base=  trapz(tbaseline,ybaseline);
        %for freq_ind = 1:size(freq_vec,2)

        colorvector = ['#0072BD';'#D95319';'#EDB120';	'#7E2F8E'];
        markervector = {'-';'--';':';'-.'};      
        kHYDbasal = kHYD;
        VIRbasal = pars(104);
for VIR_ind = 1:numel(parscan)
        for khyd_ind = 1:numel(khyd_vary)   
            pulsetime =0;
            h = khyd_vary(khyd_ind);
            kHYD = kHYDbasal * h;
            VIR = parscan(VIR_ind)*VIRbasal;
            pars(4) = kHYD;
            pars(105) = VIR;
            [t0,y_ss]= ode15s(@leung2022_SERCA,[0:1:1000],init,options,pars, kHYD,pulsetime);
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
                    [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
                    options,pars, kHYD,pulsetime);
                    y = [y;x_temp];
                    t = [t;t_temp];
                end
                recov_dur = 500;
                recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
                [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
                stimstart=recovtime(end)+0.001;
            

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
        
        t = [t0(1:end)-t_start;t];
        y = [y_ss(1:end,:);y];
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
basepars = finalpars;
pars = basepars;
pulsetime =0;
kHYD = pars(4);
VIR = pars(105);

    tstart = 1000;
    [t0,y_ss]= ode15s(@leung2022_SERCA,[0:1:tstart],init,options,pars, kHYD,pulsetime);
    initcond =  y_ss(end,:);
    [tbaseline,ybaseline]= ode15s(@leung2022_SERCA,[0:1:200],init,options,pars, kHYD,pulsetime);
    AMPK_base =  trapz(tbaseline,ybaseline);
    mTORC1_base=  trapz(tbaseline,ybaseline);
    mTORC2_base =  trapz(tbaseline,ybaseline);
    CaC_base=  trapz(tbaseline,ybaseline);
    %for freq_ind = 1:size(freq_vec,2)
    parscan = [0.5,1,2,5];
    freq_vary = [0.1,1,10,50];
    colorvector = ['#0072BD';'#D95319';'#EDB120';	'#7E2F8E'];
    markervector = {'-';'--';':';'-.'};      
    for VIR_ind = 1:numel(parscan)
        for freq_ind = 1:numel(freq_vary)    
            pulsetime =0;
            h = freq_vary(freq_ind); %vary frequency in nested for loop
            freepar = parscan(VIR_ind);
            pars(105) = freepar*VIR; % vary VIR in nested for loop
            %initializes the simulation for new parameter set
            [t0_temp,y_ss_temp]= ode15s(@leung2022_SERCA,[0:1:tstart],init,options,pars, kHYD,pulsetime); 
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
                [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 500;
            recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;
            recov_dur = 500;
            recovtime = recovtime(end)+1:recovtime(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];

            % find timepoint where t=200, etc
            [~,t200] = min(abs(t-200));
            [~,t400] = min(abs(t-400));
            [~,t600] = min(abs(t-600));
            
            
            %compute AUC's and assign to matrix
            AMPK_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,pAMPK));
            mTORC1_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,pmTORC1));
            mTORC2_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,pmTORC2));
            CaC_AUC(VIR_ind,freq_ind) = trapz(t(1:t200),y(1:t200,Ca_C));
            
            %         AMPK_AUC = AMPK_AUC ./ AMPK_base;
            %         mTORC1_AUC = mTORC1_AUC ./ mTORC1_base;
            %         mTORC2_AUC = mTORC2_AUC ./ mTORC2_base;
            %         CaC_AUC = CaC_AUC ./ CaC_base;
         
            t = [t0(1:end)-tstart;t]; % add the initialization
            y = [y_ss(1:end,:);y];
            plotfigure9
            %% Plots
        end
        
    end

subplot(3,2,5)

legend('0.5,0.1','0.5,1','0.5,10','0.5,50','1,0.1','1,1','1,10','1,50','2,0.1','2,1','2,5','2,50','5,0.1','5,1','5,10','5,50')

end
if (fig_index == 10 )
    % Frequency vs kHYD
    clear AMPK_AUC mTORC1_AUC CaC_AUC mTORC2_AUC
    %% Simulation Engine
    % initialization run
    % or load premade initialization data
    basepars = finalpars;
    pars = basepars;
    pulsetime =0;
    kHYD = pars(4);
    VIR = pars(105);
    tstart = 1000;
     freq_vary = [0.1,1,10,50];
    khyd_vary = [0.5,1,2,5];

    [~,y_ss]= ode15s(@leung2022_SERCA,[0:1:tstart],init,options,pars, kHYD,pulsetime);
    initcond =  y_ss(end,:);

    [tbaseline,ybaseline]= ode15s(@leung2022_SERCA,[0:1:200],init,options,pars, kHYD,pulsetime);

    AMPK_base =  trapz(tbaseline,ybaseline);
    mTORC1_base=  trapz(tbaseline,ybaseline);
    mTORC2_base =  trapz(tbaseline,ybaseline);
    CaC_base=  trapz(tbaseline,ybaseline);
    %for freq_ind = 1:size(freq_vec,2)

    colorvector = ['#0072BD';'#D95319';'#EDB120';	'#7E2F8E'];
    markervector = {'-';'--';':';'-.'};      
    kHYDbasal = kHYD;
    for freq_ind = 1:numel(freq_vary)
        for khyd_ind = 1:numel(khyd_vary)    
            pulsetime = 0;
            h = khyd_vary(khyd_ind);
            kHYD = kHYDbasal * h;
            frequency = freq_vary(freq_ind);
            [t0,y_ss]= ode15s(@leung2022_SERCA,[0:1:tstart],init,options,pars, kHYD,pulsetime);
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
                [t_temp,x_temp]= ode15s(@leung2022_SERCA,[startpoint endpoint],initcond,...
                options,pars, kHYD,pulsetime);
                y = [y;x_temp];
                t = [t;t_temp];
            end
            recov_dur = 500;
            recovtime = timevec(end)+1:1:timevec(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];
            stimstart=recovtime(end)+0.001;
            
            recov_dur = 500;
            recovtime = recovtime(end)+1:recovtime(end) +recov_dur;
            [t_temp,x_temp]= ode15s(@leung2022_SERCA,[recovtime],y(end,:),...
            options,pars, kHYD,pulsetime);
            y = [y;x_temp];
            t = [t;t_temp];

            hold on
            [~,t200] = min(abs(t-200));
            [~,t400] = min(abs(t-400));
            [~,t600] = min(abs(t-600));
            
            
            
            AMPK_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pAMPK));
            mTORC1_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pmTORC1));
            mTORC2_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,pmTORC2));
            CaC_AUC(freq_ind,khyd_ind) = trapz(t(1:t200),y(1:t200,Ca_C));
            

            t = [t0(1:end)-tstart;t];
            y = [y_ss(1:end,:);y];
            plotfigure10
            %% Plots
        end
    end

subplot(3,2,5)
     freq_vary = [0.1,1,10,50];
    khyd_vary = [0.5,1,2,5];
legend('0.1,0.5','0.1,1','0.1,2','0.1,5','1,0.5','1,1','1,2','1,5','10,0.1','10,1','10,2','10,5','50,0.1','50,1','50,2','50,5')
end


 %       legend('1','2','3','4','5','6','7')
end