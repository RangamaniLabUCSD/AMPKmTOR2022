
function [f] = leung2022_SERCA(t,x,pars,kHYD,pulsetime)
        %==================================================================
%% Species Definition
    %ATP module
        ATP       = x(1);
        ADP       = x(2);
        AMP       = x(3);
        PCr       = x(4);
        Pi        = x(5);
        AMPK      = x(6);
        pAMPK     = x(7);
       Act_AMPK  = x(8);
        Act_pAMPK = x(9);
        Act       = x(10);
        AICAR     = x(11);
    % Calcium module
        Glut      =  x(12);
        Ca_ER     =  x(13);
        Ca_C      =  x(14); %
        w         =  x(15);%15
        Ri        =  x(16);  
        R2        =  x(17);
        DIM       =  x(18);
        DAG       =  x(19);
        DIMp      =  x(20);%20
        NMDA_C0   =  x(21);
        NMDA_C1   =  x(22);
        NMDA_C2   =  x(23);
        NMDA_D    =  x(24);
        NMDA_O    =  x(25);%25
        PKC       =  x(26);
        IP3       =  x(27);
        B         =  x(28);
        BCa       =  x(29);
        BuffER    =  x(30);%30
        AMPA_U    =  x(31); 
        AMPA_M    =  x(32);
        AMPA_C    =  x(33);
        AMPA_O    =  x(34); 
        AMPA_D1   =  x(35);%35
        AMPA_D2   =  x(36);
        AMPA_D3   =  x(37);
        
        V         = x(38);
        CaM       =  x(39);
        CaCaM     =  x(40);%40
        CaMKK2    =  x(41);
        CaMKK2_act=  x(42);
        %% MTOR AMPK
    IR            = 1000 * x(43);
    pIR           = 1000 * x(44);
    IRS           = 1000 * x(45);
    pIRS          = 1000 * x(46);
    iIRS          = 1000 * x(47);
    AKT           = 1000 * x(48);
    pAKT          = 1000 * x(49);
    mTORC1        = 1000 * x(50);
    pmTORC1       = 1000 * x(51);
    mTORC2        = 1000 * x(52);
    pmTORC2       = 1000 * x(53);
    mTORC1_DEPTOR = 1000 * x(54);
    mTORC2_DEPTOR = 1000 * x(55);
    DEPTOR        = 1000 * x(56);
    pDEPTOR       = 1000 * x(57);
    SIRT1         = 1000 * x(58);
    ULK1          = 1000 * x(59);
    pULK1         = 1000 * x(60);




        
%% Parameter Definition
    %ATP/AMPK module
        KADP        = pars(1);
        VmaxOP      = pars(2);
        nH          = pars(3);

        krest      = pars(4);
        if(kHYD>0)
            krest = kHYD;
        end
        kstim      = pars(5);
        kpost      = pars(6);
        vAK         = pars(7);
        kmt         = pars(8);
        kmd         = pars(9);
        kmm         = pars(10);
        keqadk      = pars(11);
       
        vmax20      = pars(12);
        vmax21      = pars(13);
        km20        = pars(14);
        km21        = pars(15);
        k12f        = pars(16);
        k12r        = pars(17);
        k13f        = pars(18);
        k13r        = pars(19);
        kaicar      = pars(20);
        kAct        = pars(21);
        kcatAMPK    = pars(22);
        VforCK      = pars(23);
        Kb          = pars(24);
        Kia         = pars(25);
        Kib         = pars(26);
        Kiq         = pars(27);
        Kp          = pars(28);
        KeqCK       = pars(29);
        TCr         = pars(30);
        
 %% Calcium module       
        %IP3 and mgluR
        k1          = pars(31);% 30;%1/s
        b           = pars(32);% 0.001;
        Ki          = pars(33);% 5; %uM
        Ka          = pars(34);% 3; %uM
        IRa         =  (1-Ri)*(IP3^2/(Ki^2+IP3^2))*(Ca_C^3/(Ka^3+Ca_C^3));
        k_          = pars(36);% 0.02; % 1/s
        k           = pars(37);% 20;   % % 1/uM^4s
        VDAG        = pars(38);% 0.0325 ; %1/uM
        kb          = pars(39);% 0.1; %1/uM^2s
        ku          = pars(40);% 2; %
        Vp          = pars(41);% 0.05; %UM/s
        Kmp         = pars(42);% 5e-4; %uM
        Vpkc        = pars(43);% 0.2; % uM/s
        Kmpkc       = pars(44);% 5e-4; %uM
        kapkc       = pars(45);% 0.2; %1/s
        kdpkc       = pars(46);% 20; % 1/s
        kplcI       = pars(47);% 3; %1/s
        kplcD       = pars(48);% 3; %1/s
        KmDAG       = pars(49);% 6; %uM
        
        
        %RYR
        VRYR        = pars(50);% 0.0050; %1/s
        Kar         = pars(51);% 0.0192;%1/s
        Kbr         = pars(52);% 0.2573;%1/s
        Kcr         = pars(53);% 0.0571;%1/s
        Kdr         = pars(54);% 0.1;%1/s
        
        %SERCA
        VSERCA      = pars(55);% 120;
        Kp          = pars(56);% 6;
        Ke_SERCA = 2000;
        %NMDA
        Rb          = pars(57);% 5; % 1/uM.s
        Ru          = pars(58);% 46.5; % /s
        Rd          = pars(59);% 8.4; % 1/s
        Rr          = pars(60);% 6.8; % 1/s
        Ro          = pars(61);% 46.5; % 1/s
        Rc          = pars(62);% 73.8; % 1/s
        tauess      = pars(63);% 0.05; % s
        tauesf      = pars(64);% 0.005; % s
        taubss      = pars(65);% 0.025; % s
        taubsf      = pars(66);% 0.003; % s
        tdelaybp    = pars(67);% 0.002; %s
        s_term      = pars(68);% 10; % mV
        V_reversal  = -pars(69);% (-)65; % mV
        N_NMDA      = pars(70);% 1; 
        BPAPmax     = pars(71);% 40; %mV
        G_NMDA      = -pars(72);% (-)65.6e6;
        k_ext       = pars(73);%0.25*1e19/(2 * 1.602*6.022e23);%picoAmpere to Flux (micro mol/Ls)  
        Kmext       = pars(74);% 0.1;
        kdeg        = pars(75);% 5;
        kpmleak     = pars(76);% 0.5;
        k_buffcyt   = pars(77);%   0.1;
        k_buffcyt_r = pars(78);%   1 ;
        k_ERBUFF    = pars(79);%  50;
        k_ERBUFF_   = pars(80);%  3 ;
        kerleak     = pars(81);% 0.1; 
 %AMPA
        G_AMPA      = pars(82);% 15;
        kAMPA_1f    = pars(83);% 1.8 ; % 1/uMs
        kAMPA_2f    = pars(84);% 10;
        kAMPA_3f    = pars(85);% 1.6e4;
        kAMPA_4f    = pars(86);% 7e2;
        kAMPA_5f    = pars(87);% 1e2;
        kAMPA_6f    = pars(88);% 3e2;
        kAMPA_7f    = pars(89);% 10;
        kAMPA_8f    = pars(90);% 1.6e4;
        
        kAMPA_1r    = pars(91);% 2.4e3;
        kAMPA_2r    = pars(92);% 1e4;
        kAMPA_3r    = pars(93);% 5e3;
        kAMPA_4r    = pars(94);% 1.5e2;
        kAMPA_5r    = pars(95);% 2.1;
        kAMPA_6r    = pars(96);% 15;
        kAMPA_7r    = pars(97);% 1e3;
        kAMPA_8r    = pars(98);% 1.2e4;

        k_cam_f     = pars(99);% 0.001;
        k_cam_r     = pars(100);% 01.0;
        Vm_kk2      = pars(101);%0.01;
        Km_kk2      = pars(102);%0.01;
        Ka_kk2      = pars(103);%0.005;
        kf_CaMDisc  = pars(104);% 1;

%% mTOR AMPK
       parstart = 104;
      
     V_IR                        = pars(parstart +1);         % Rate of activation of IR
    Km_IR                       = pars(parstart +2);         % MM constant for the activation of IR
    V_pIR                       = pars(parstart +3);         % Rate of deactivation of IR
    Km_pIR                      = pars(parstart +4);         % MM constant for the deactivation of IR
    K_IRS_by_pIR                = pars(parstart +5);         % Rate of activation of IRS via pIR
    Km_IRS_by_pIR               = pars(parstart +6);         % MM constant for the activation of IRS via pIR
    V_pIRS                      = pars(parstart +7);         % Rate of deactivation of IRS
    Km_pIRS                     = pars(parstart +8);         % MM constant for the deactivation of IRS
    K_AKT_by_pIRS               = pars(parstart +9);         % Rate of activation of AKT via pIRS
    Km_AKT_by_pIRS              = pars(parstart +10);        % MM constant for the activation of AKT via pIRS
    K_AKT_by_pmTORC2            = pars(parstart +11);        % Rate of activation of AKT via pmTORC2
    Km_AKT_by_pmTORC2           = pars(parstart +12);        % MM constant for the activation of AKT via pmTORC2
    V_pAKT                      = pars(parstart +13);        % Rate of deactivation of AKT
    Km_pAKT                     = pars(parstart +14);        % MM constant for the deactivation of AKT
    K_mTORC1_by_pAKT            = pars(parstart +15);        % Rate of activation of mTORC1 via pAKT
    Km_mTORC1_by_pAKT           = pars(parstart +16);        % MM constant for the activation of mTORC1 via pAKT
    K_pmTORC1                   = pars(parstart +17);        % Rate of background deactivation of mTORC1
    K_pmTORC1_by_pAMPK          = pars(parstart +18);        % Rate of deactivation of mTORC1 via pAMPK
    Km_pmTORC1_by_pAMPK         = pars(parstart +19);        % MM constant for the deactivation of mTORC1
    K_pmTORC1_by_pULK1          = pars(parstart +20);        % Rate of deactivation of mTORC1 via pULK1
    Km_pmTORC1_by_pULK1         = pars(parstart +21);        % MM constant for the deactivation of mTORC1 via pULK1
    K_mTORC2_by_pIRS            = pars(parstart +22);        % Rate of activation of mTORC2 via pIRS
    Km_mTORC2_by_pIRS           = pars(parstart +23);        % MM constant for the activation of mTORC2 via pIRS
    K_mTORC2_by_pAMPK           = pars(parstart +24);        % Rate of activation of mTORC2 via pAMPK
    Km_mTORC2_by_pAMPK          = pars(parstart +25);        % MM constant for the activation of mTORC2 via pAMPK
    V_pmTORC2                   = pars(parstart +26);        % Rate of deactivation of mTORC2
    Km_pmTORC2                  = pars(parstart +27);        % MM constant for the deactivation of mTORC2
    K_DEPTOR_by_pmTORC1         = pars(parstart +28);        % Rate of activation of DEPTOR via pmTORC1
    Km_DEPTOR_by_pmTORC1        = pars(parstart +29);        % MM constant for the activation of DEPTOR via pmTORC1
    K_DEPTOR_by_pmTORC2         = pars(parstart +30);        % Rate of activation of DEPTOR via pmTORC2
    Km_DEPTOR_by_pmTORC2        = pars(parstart +31);        % MM constant for the activation of DEPTOR via pmTORC2
    V_pDEPTOR                   = pars(parstart +32);        % Rate of deactivation of DEPTOR
    Km_pDEPTOR                  = pars(parstart +33);        % MM constant for the deactivation of DEPTOR
    K_mTORC1_DEPTOR_form        = pars(parstart +34);        % Rate of formation of the mTORC1-DEPTOR complex
    K_mTORC1_DEPTOR_diss        = pars(parstart +35);        % Rate of dissociation of the mTORC1-DEPTOR complex
    K_mTORC2_DEPTOR_form        = pars(parstart +36);        % Rate of formation of the mTORC2-DEPTOR complex
    K_mTORC2_DEPTOR_diss        = pars(parstart +37);        % Rate of dissociation of the mTORC2-DEPTOR complex
    K_IRS_to_iIRS               = pars(parstart +38);        % Rate of inactivation of IRS
    Km_IRS_to_iIRS              = pars(parstart +39);        % MM constant for the inactivation of IRS
    V_iIRS                      = pars(parstart +40);        % Rate of activation of IRS from iIRS
    Km_iIRS                     = pars(parstart +41);        % MM constant for the activation of IRS from iIRS
    K_AMPK                      = pars(parstart +42);        % Rate of background activation of AMPK
    K_AMPK_by_SIRT1             = pars(parstart +43);        % Rate of activation of AMPK via SIRT1
    Km_AMPK                     = pars(parstart +44);        % MM constant for the activation of AMPK
    K_pAMPK                     = pars(parstart +45);        % Rate of background deactivation of AMPK
    K_pAMPK_by_pULK1            = pars(parstart +46);        % Rate of deactivation of AMPK via pULK1
    K_pAMPK_by_pmTORC1          = pars(parstart +47);        % Rate of deactivation of AMPK via pmTORC1
    Km_pAMPK                    = pars(parstart +48);        % MM constant for the deactivation of AMPK
    K_SIRT1                     = pars(parstart +49);        % Rate of background activation of SIRT1
    K_SIRT1_by_pAMPK            = pars(parstart +50);        % Rate of activation of SIRT1 via pAMPK
    Km_SIRT1                    = pars(parstart +51);        % MM constant for the activation of SIRT1
    K_SIRT1_diss                = pars(parstart +52);        % Rate of dissociation of SIRT1
    K_ULK1                      = pars(parstart +53);        % Rate of background activation of ULK1
    K_ULK1_by_pAMPK             = pars(parstart +54);        % Rate of activation of ULK1 via pAMPK
    Km_ULK1                     = pars(parstart +55);        % MM constant for the activation of ULK1
    K_pULK1                     = pars(parstart +56);        % Rate of background deactivation of pULK1
    K_pULK1_by_pmTORC1          = pars(parstart +57);        % Rate of deactivation of pULK1 via pmTORC1
    Km_pULK1                    = pars(parstart +58);        % MM constant for the deactivation of ULK1

    krampk_bal = pars(parstart+59);
    km_CA_ATP                   = 500;

  %% Flux Definition : Mitochondrial and Cytosolic metabolism     
        
        
        %modified to have krest be the varying consumption rate
%% test
   %krest = multkrest*krest;
   %k7f = k7f/50;
   %k10f = k10f/50;
  %  KpmTORC1AMPK= multkrest*KpmTORC1AMPK;
%krampk_bal = multkrest;
%IRS = multkrest;

%%

        r1  = krest*ATP;
   
        r2_num=VmaxOP*(ADP/KADP)^nH;%^*(Ca_C^1.4/(Ca_C^1.4 + km_CA_ATP));
        r2_den=(1+(ADP/KADP)^nH);
        r2  = r2_num/r2_den;
        
        r3_num=VforCK * ADP * PCr / (Kia*Kb);
        r3_den=(1+ADP/Kia+PCr/Kib+ATP/Kiq+(ADP*PCr/(Kia*Kb))+((TCr - PCr)* ATP /(Kiq*Kp)));
        r3  =  r3_num/r3_den;
        
        vCKr = VforCK * Kiq * Kp/ (KeqCK * Kia *Kb);
        r4_num=vCKr * ATP *(TCr - PCr)/(Kiq*Kp);
        r4_den=(1+(ADP/Kia)+(PCr/Kib)+(ATP/Kiq)+(ADP*PCr/(Kia*Kb))+((TCr - PCr)* ATP /(Kiq*Kp)));
        r4  =  r4_num/r4_den;
        
        CK = r4-r3;
        r5f_num=(vAK*ATP*(AMP)/(kmt*kmm));
        r5f_den=(1+(ATP/kmt)+(AMP/kmm)+(ATP*AMP/(kmt*kmm))+(2*ADP/kmd)+(ADP^2/kmd^2));
        r5f = r5f_num/r5f_den;
        
        VrevAK=(vAK*kmd^2/(keqadk*kmt*kmm));
        
        r5r_num=VrevAK*ADP^2/kmd^2;
        r5r_den=(1+(ATP/kmt)+(AMP/kmm)+(ATP*AMP/(kmt*kmm))+(2*ADP/kmd)+(ADP^2/kmd^2));
        r5r = r5r_num/r5r_den;
        AK =  r5f-r5r;


        rc = 1300* (AMP * AMPK - krampk_bal/1300 *ATP * pAMPK);



        r20 = vmax20*Act_AMPK/(km20 + Act_AMPK);
        r21 = vmax21*Act_pAMPK/(km21 + Act_pAMPK ); 
        r22 = k12f*AMPK*Act - k12r*Act_AMPK;
        r23 = k13f*pAMPK*Act -k13r * Act_pAMPK;
        r24 = kaicar * AICAR;
        r25 = kAct * Act;


  %% Calcium Module
        %IP3
        J_Ri = k*Ca_C^4*(1-Ri)/(1+(Ca_C/Ka)^3)-k_*Ri;
      
        J_IP3R = k1*(b+IRa)*(Ca_ER);   
        % RYR
        %
        Po = w*(1+(Ca_C^3/Kbr))/((Kar/Ca_C^4)+1+(Ca_C^3/Kbr)); 
        winf = ((Kar/Ca_C^4)+1+(Ca_C^3/Kbr))/((1/Kcr)+(Kar/Ca_C^4)+1+(Ca_C^3/Kbr));
        tau = winf/Kdr;
        J_w = (winf-w)/tau;
        J_RYR =  VRYR * Po*(Ca_ER); 

        JPMCA = k_ext * (Ca_C)^2/((Ca_C)^2 + Kmext)*(ATP)/(Ke_SERCA+ATP);

        J_PMleak =1* kpmleak * (Ca_C - 0.1);
       
        
        %SERCA
        %kxserca = (k2serca*(Ca_C^2)+ k4rserca * K1serca)/(K1serca + Ca_C^2);
        %kyserca = (k2rserca * K3serca*(Ca_ER^2) +k4serca)/(gammaserca *(1+K3serca*Ca_ER^2));
        %J_SERCA = 
        J_SERCA = VSERCA*(Ca_C^2/(Kp^2+Ca_C^2))*(ATP)/(Ke_SERCA+ATP);
        %

        J_BuffCa = k_buffcyt * Ca_C * B -k_buffcyt_r * BCa;

        J_ER_Buff = k_ERBUFF * Ca_ER - k_ERBUFF_ *BuffER;
       
        J_ER_leak = 1*kerleak*(Ca_ER - 600);
        %NMDA
        JN1   = Rb*NMDA_C0*Glut - Ru*NMDA_C1;
        JN2   = Rb*NMDA_C1*Glut- Ru*NMDA_C2;
        JN3   = Rd*NMDA_C2 - Rr*NMDA_D;
        JN4   = Ro*NMDA_C2 - Rc*NMDA_O;
        BPAP =BPAPmax*(0.75*exp(-(t-pulsetime -tdelaybp)/taubsf) +0.25*exp(-(t-pulsetime)/taubss)  )/2;
        EPSP =s_term *(0.5*exp(-(t-pulsetime)/tauesf) +  0.5*exp(-(t-pulsetime)/tauess));
        BV = 1/(1 + (exp(-0.092*V ) * (1/3.57)));
                
 


        
        
        JA1 = kAMPA_1f *AMPA_U * Glut -kAMPA_1r *AMPA_M;
        JA2 = kAMPA_2f *AMPA_M * Glut-kAMPA_2r *AMPA_C;
        JA3 = kAMPA_3f *AMPA_C -kAMPA_3r *AMPA_O;
        JA4 = kAMPA_4f *AMPA_M -kAMPA_4r *AMPA_D1;
        JA5 = kAMPA_5f *AMPA_C -kAMPA_5r *AMPA_D2;
        JA6 = kAMPA_6f *AMPA_O -kAMPA_6r *AMPA_D3;
        JA7 = kAMPA_7f *AMPA_D1* Glut -kAMPA_7r *AMPA_D2;
        JA8 = kAMPA_8f *AMPA_D2 -kAMPA_8r *AMPA_D3;
        AMPA_EPSP = AMPA_O * G_AMPA *BV* (V);

        V = -65 + BPAP + EPSP + AMPA_EPSP;
        AMPA_EPSP = AMPA_O * G_AMPA *BV* (V);
         BPAP =0.5*BPAPmax*(0.75*exp(-(t-pulsetime -tdelaybp)/taubsf) +0.25*exp(-(t-pulsetime)/taubss)  );
        EPSP =s_term *(0.5*exp(-(t-pulsetime)/tauesf) +  0.5*exp(-(t-pulsetime)/tauess));
        BV = 1/(1 + (exp(-0.092*V ) * (1/3.57)));
                J_NMDA = 0.1*G_NMDA*NMDA_O*BV*((V+V_reversal))*N_NMDA*1e19/ (2 * 1.602*6.022e23);      

        J_PM = 5*(J_NMDA  - JPMCA - J_PMleak) ;

        %mgluR  
        J_R2    =  -kb*R2*Glut^2+ku*DIM;
        J_DIM   =  kb*R2*Glut^2-ku*DIM+Vp*DIMp/(Kmp+DIMp)-Vpkc*PKC*DIM/(Kmpkc+DIM);
        J_DIMp  =  -Vp*DIMp/(Kmp+DIMp)+Vpkc*PKC*DIM/(Kmpkc+DIM);
        J_PKC   =  kapkc*DAG/(DAG+KmDAG)*(1-PKC)-kdpkc*PKC;
        J_IP    =  kplcI*DIM;
        J_DAG   =  kplcD*DIM-VDAG*DAG/(DAG+KmDAG);
        
        
        J_CaMBind = k_cam_f * Ca_C * CaM;
        J_CaMDisc = k_cam_r * CaCaM;
        J_CAMKK2_Act = Vm_kk2 * (1 + CaCaM/Ka_kk2)*CaMKK2/(Km_kk2 + CaMKK2);
        J_CAMKK2_Deac = kf_CaMDisc * CaMKK2_act;
        
        J_ATP_Ca = J_SERCA + JPMCA;



SIRT1_total = 375;
%pre convert for JM fluxes, which are in units of mM
AMPK = (Act_AMPK + AMPK)*1000;
pAMPK = (Act_pAMPK + pAMPK)*1000;


       % J_ATP_Ca = 0;

%% parbox
%    KpmTORC1 = KpmTORC1 * 1;
   %  KmpmTORC1A   = KmpmTORC1A * parmult;
%  




    %% Fluxes for mTOR
        % format FluxName = mm(Vmax,Km, Substrate)
        % Output
    JM1  = (V_IR*IR)/(Km_IR+IR); % IR...pIR
    JM2  = (V_pIR*pIR)/(Km_pIR+pIR); % pIR...IR
    JM3  = (K_IRS_by_pIR*pIR*IRS)/(Km_IRS_by_pIR+IRS); % IRS...pIRS
    JM4  = (V_pIRS*pIRS)/(Km_pIRS+pIRS); % pIRS...IRS
    JM5  = (K_AKT_by_pIRS*pIRS*AKT)/(Km_AKT_by_pIRS+AKT)+(K_AKT_by_pmTORC2*pmTORC2*AKT)/(Km_AKT_by_pmTORC2+AKT); % AKT...pAKT
    JM6  = (V_pAKT*pAKT)/(Km_pAKT+pAKT); % pAKT...AKT
    JM7  = (K_mTORC1_by_pAKT*pAKT*mTORC1)/(Km_mTORC1_by_pAKT+mTORC1); % mTORC1...pmTORC1
    JM8  = (K_pmTORC1 + K_pmTORC1_by_pAMPK*pAMPK)*pmTORC1/(Km_pmTORC1_by_pAMPK+pmTORC1)+(K_pmTORC1_by_pULK1*pULK1)*pmTORC1/(Km_pmTORC1_by_pULK1+pmTORC1); % pmTORC1....mTORC1
    JM9  = (K_mTORC2_by_pIRS*pIRS*mTORC2)/(Km_mTORC2_by_pIRS+mTORC2)+(K_mTORC2_by_pAMPK*pAMPK*mTORC2)/(Km_mTORC2_by_pAMPK+mTORC2); % mTORC2...pmTORC2
    JM10 = (V_pmTORC2*pmTORC2)/(Km_pmTORC2+pmTORC2); % pmTORC2...mTORC2
    JM11 = (K_DEPTOR_by_pmTORC1*pmTORC1*DEPTOR)/(Km_DEPTOR_by_pmTORC1+DEPTOR)+(K_DEPTOR_by_pmTORC2*pmTORC2*DEPTOR)/(Km_DEPTOR_by_pmTORC2+DEPTOR); % DEPTOR...pDEPTOR
    JM12 = (V_pDEPTOR*pDEPTOR)/(Km_pDEPTOR+pDEPTOR); % pDEPTOR...DEPTOR
    JM13 = (K_mTORC1_DEPTOR_form*mTORC1*DEPTOR)-(K_mTORC1_DEPTOR_diss*mTORC1_DEPTOR); % mTORC1+DEPTOR...mTORC1_DEPTOR
    JM14 = (K_mTORC2_DEPTOR_form*mTORC2*DEPTOR)-(K_mTORC2_DEPTOR_diss*mTORC2_DEPTOR); % mTORC2+DEPTOR...mTORC2_DEPTOR 
    JM15 = (K_IRS_to_iIRS*pmTORC1*IRS)/(Km_IRS_to_iIRS+IRS); % IRS...iIRS
    JM16 = (V_iIRS*iIRS)/(Km_iIRS+iIRS); % iIRS...IRS
    JM17 = (K_AMPK+K_AMPK_by_SIRT1 * SIRT1)*AMPK/(Km_AMPK + AMPK); % AMPK...pAMPK
    JM18 = (K_pAMPK+K_pAMPK_by_pULK1 * pULK1 + K_pAMPK_by_pmTORC1 * pmTORC1)*pAMPK/(Km_pAMPK + pAMPK); % pAMPK...AMPK
    JM19 = (K_SIRT1+K_SIRT1_by_pAMPK * pAMPK)*(SIRT1_total - SIRT1) /(Km_SIRT1 + SIRT1_total - SIRT1) - K_SIRT1_diss * SIRT1; % pAMPK...SIRT
    JM20 = (K_ULK1+K_ULK1_by_pAMPK * pAMPK)*ULK1/(Km_ULK1 + ULK1); % ULK1...pULK1
    JM21 = (K_pULK1+K_pULK1_by_pmTORC1 * pmTORC1)*pULK1/(Km_pULK1 + pULK1); % pULK1...ULK1
AMPK = (AMPK-Act_AMPK*1000)/1000;
pAMPK = (pAMPK-Act_pAMPK*1000)/1000; % reconvert back to normal (this doesn't have an explicit function right now, but is useful for later

    %% ODE Definition


        

  %%
    



        f(1)  = -r1 + r2 -AK -CK  -2* J_ATP_Ca +rc/1000;%ATP
        f(2)  = r1-r2 + 2*AK +CK +2*J_ATP_Ca ;%ADP
        f(3)  = -AK  -r20 - rc/1000;%AMP
        f(4)  = CK;
        f(5)  = r1 -r2;%Pi
        f(5)  =0; % clamp Pi, it is not used anyway
        f(6) = -rc + JM18-JM17;%AMPK
        f(7) = rc - JM18 + JM17;%pAMPK
        f(8) = -r20-r21+r22;%Act_AMPK
        f(9) = r20-r21+r23;%Act_pAMPK
        f(10) = -r22-r23+r24-r25;%Act
        f(11) = 0;%AICAR

        f(8:11) = 0;
        f(12) =  -200 * Glut; % Glutamate
        f(13) =  4*(J_SERCA - J_IP3R - J_RYR-J_ER_leak) - J_ER_Buff; % Ca_ER
        %xdot(2) = 0;
        f(14) =  -1/4*(J_SERCA - J_IP3R - J_RYR+J_ER_leak) + J_PM - J_BuffCa -J_CaMBind; % CA_C
        f(15) =  J_w; % w
        f(16) =  J_Ri; % Ri
        f(17) =  J_R2; % R2
        f(18) =  J_DIM; % DIM
        f(19) =  J_DAG; %DAG
        f(20) =  J_DIMp; %dimp
        f(21) = -JN1 ; %NMDA C0
        f(22) =  +JN1 - JN2; %NMDA C1
        f(23) =  +JN2 -JN3 -JN4; %NMDA C2
        f(24) =  JN3; %NMDA D
        f(25) =  JN4; %NMDA O
        f(26) =  J_PKC; %PKC
        f(27) = J_IP - kdeg*(IP3 - 0.01);
        f(28) = -J_BuffCa;
        f(29) = +J_BuffCa;
        f(30) = J_ER_Buff;
        
        f(31) = -JA1 ;%U
        f(32) = +JA1 - JA2-JA4 ;%M
        f(33) = JA2 - JA3 - JA5 ;%C
        f(34) =  JA3 - JA6;%O
        f(35) =  JA4 - JA7;%D1
        f(36) =  JA5 + JA7 - JA8;%D2
        f(37) =  JA6 + JA8;%D3
        f(38) =  0; %V
        f(39) = -J_CaMBind + J_CaMDisc; %CaM
        f(40) = J_CaMBind - J_CaMDisc;%CaCaM
        f(41) =  -J_CAMKK2_Act + J_CAMKK2_Deac; %CAMKK2
        f(42) = J_CAMKK2_Act -J_CAMKK2_Deac; %CAMKK2_act
        f(43) = JM2-JM1;
        f(44) = JM1-JM2;
        f(45) = JM4+JM16-JM3-JM15;
        f(46) = JM3-JM4;
        f(47) = JM15-JM16;
        f(48) = JM6-JM5;
        f(49) = JM5-JM6;
        f(50) = JM8-JM7-JM13;
        f(51) = JM7-JM8;
        f(52) = JM10-JM9-JM14;
        f(53) = JM9-JM10;
        f(54) = JM13;
        f(55) = JM14;
        f(56) = JM12-JM11-JM13-JM14;
        f(57) = JM11-JM12;
        f(58) = JM19; 
        f(59) = JM21-JM20; 
        f(60) = JM20-JM21; 
        f(43:60) = f(43:60) ./1000;


        f=f';
        
end

function V = mm(Vmax,km,Substrate)
    V = Vmax*(Substrate)/(km + Substrate);
end