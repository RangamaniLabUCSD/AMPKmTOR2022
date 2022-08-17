function [coinit] = leung_icsource()
coinit = zeros(1,55);

coinit(01) = 2.6;
coinit(02) = 0.3;
coinit(03) = 0.045;
coinit(04) = 22.1;
coinit(05) = 3.183;
coinit = coinit.*1000; % convert from mm to muM

coinit(06) = .200279483201148 ;% AMPK
coinit(07) =  .05072051679885231 ;% pAMPK
coinit(08) = 0;
coinit(09) = 0;
coinit(10) = 0;
coinit(11) = 0;

coinit(12)  = 0;% Glutamate
coinit(13)  = 600;
coinit(14)  = 0.1;
coinit(15)  = 0.963;
coinit(16)  = 0;
coinit(17)  = 0.03;
coinit(18)  = 0.014;
coinit(19)  = 0.025;
coinit(20)  = (0.075-(0.1*0.002)^(1/2)-2*0.002-2*0.014)/2;

coinit(21)  = 1;
coinit(22)  = 0;
coinit(23)  = 0;
coinit(24)  = 0;
coinit(25)  = 0;
coinit(26)  = 0.2;
coinit(27)  = 0.1 ;
coinit(28)  =  20  ;
coinit(29)  =  5  ;
coinit(30)  =  1000  ;

coinit(31)  = 0.1 ; %AMPA_U
coinit(32)   =  0  ; %AMPA_M
coinit(33)   =  0  ; %AMPA_C
coinit(34)   =  0  ; %AMPA_O
coinit(35)   =  0  ; %AMPA_D1
coinit(36)   =  0  ; %AMPA_D2
coinit(37)   =  0  ;%AMPA_D3
coinit(38)   =  -65 ; %V
coinit(39)        = 10    ;% CaM
coinit(40)        =  0   ;% Ca_CaM

coinit(41)        = 10    ;%CaMKK2
coinit(42)        =  0   ;%CaMKK2_act

coinit(43)         =  45.2679        ;% IR
coinit(44)         =  4.7321        ;%pIR
coinit(45)         =  17.9913        ;%IRS
coinit(46)         =  10.4134        ;%pIRS
coinit(47)         =  71.5953        ;%iIRS

coinit(48)         =  28.9130        ;%AKT
coinit(49)         =  500      ;% pAKt
coinit(50)         =  500        ;%mTORC1
coinit(51)         =  0.9       ;%pmTORC1
coinit(52)         = 0.9      ;%mTORC2


%% Set Deptor to 0
%coinit(49:52) = 0;



coinit(53)         =  5.0290        ;% pmTORC2
coinit(54)         =  152.9209        ;%mtorC1_DEPTOR
coinit(55)         =  183.6542        ;% mTORC2_Deptor
coinit(56)         =  13.1621        ;% deptor
coinit(57)         =  0.2628        ;% pdeptor

coinit(58)         =  37.0602        ;%SIRT
coinit(59)         =  199.6328        ;% ULK
coinit(60)         =  50.3672        ;%\pULK
% y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;0;250;0]; %Initial condition
  y0 = [50;1;100;1;1;...
      100;1;250;10;200;...
      10;1;1;350;1;...
      1;250;1]; %Initial condition

coinit(43:60) = y0;
coinit(43:60) = coinit(43:60).* 1e-3; % conversion from nM to uM

% y0 = [50;0;100;0;0;100;0;250;0;200;0;0;0;350;0;0;250;0]; %Initial condition
% y0 = [45.2678886688302,4.73211133116980,17.9912537197180,10.4134176890683,71.5953285901381,...
% 28.9130204869018,71.0869795130983,87.9364799719260,9.14262282859870,11.3168282120687,5.02897769163197,152.920897222202,...
% 183.654194080837,13.1620960677460,0.262812643184197,37.0602139519015,199.632839162159,50.3671608378409];

%  %Define each network componenet to the corresponding vector element
%     IR            = x(1);
%     pIR           = x(2);
%     IRS           = x(3);
%     pIRS          = x(4);
%     iIRS          = x(5);
%     AKT           = x(6);
%     pAKT          = x(7);
%     mTORC1        = x(8);
%     pmTORC1       = x(9);
%     mTORC2        = x(10);
%     pmTORC2       = x(11);
%     mTORC1_DEPTOR = x(12);
%     mTORC2_DEPTOR = x(13);
%     DEPTOR        = x(14);
%     pDEPTOR       = x(15);
%     AMPK          = x(16);
%     pAMPK         = x(17);
%     SIRT1         = x(18);
%     ULK1          = x(19);
%     pULK1         = x(20);
end

