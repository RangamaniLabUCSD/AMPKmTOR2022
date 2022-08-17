function [copars] = leung_pars_source()
copars = zeros(1,104);
%% Metabolic model - Coccimiglio
copars(01) = 500;   % KADP
copars(02) = 0.5*1000;      % VmaxOP
copars(03) = 2.568;         % nH
copars(04) = 2.6e-2;        % krest
copars(05) = 5.0e-2;        % kstim
copars(06) = 2.6e-2;        % kpost
copars(07) = 14.66*1000;    % vAK
copars(08) = 0.27*1000;     % kmt
copars(09) = 0.35*1000;     % kmd


copars(10) = 0.32*1000;     % kmm
copars(11) = 0.744;         % keqadk

copars(12) = 2e-2*1000;     % vmax20
copars(13) = 1e-4*1000;     % vmax21
copars(14) = 1.4*1000;      % km20
copars(15) = 6.7e-2*1000;   % km21
copars(16) = 1/1000;        % k12f
copars(17) = 1.8e-2;        % k12r
copars(18) = 1/1000;        % k13f
copars(19) = 1.8e-2;        % k13r
copars(20) = 4e-4;          % kaicar
copars(21) = 3e-3;          % kAct

copars(22) = 34.56;         % kcatAMPK
copars(23) = 1e2 *1000 ;    % VCK
copars(24) = 1.11  *1000;   % Kb
copars(25) = 0.135*1000  ;  % Kia
copars(26) = 3.9 *1000 ;    % Kib
copars(27) = 3.5 *1000 ;    % Kiq
copars(28) =  3.8*1000 ;    % Kp
copars(29) = 1.77*10^(9-7); % KeqCK
copars(30) =  39*1000 ;     % TCr

%% Calcium Module Parameters
copars(31) =        30;     %1/s  k1

copars(32) =       0.001;   %     b
copars(33) =       5;       %uM   Ki
copars(34) =        3;      %u M  Ka
copars(35) =       0;            % IRA (specified in model)
copars(36) =        0.02;   % 1/s
copars(37) =        20;     % % 1/uM^4s
copars(38) =       0.0325 ; %1/uM
copars(39) =       0.1;     %1/uM^2s
copars(40) =       2;       %
copars(41) =       0.05;    %UM/s

copars(42) =       5e-4;    %uM
copars(43) =       0.2;     % uM/s
copars(44) =      5e-4;     %uM
copars(45) =      0.2;      %1/s
copars(46) =        20;     % 1/s
copars(47) =         3;     %1/s
copars(48) =    3;          %1/s
copars(49) =      6;        %uM
copars(50) =  0.0050;       %1/s
copars(51) =    0.0192;     %1/s

copars(52) =      0.2573;   %1/s
copars(53) =      0.0571;   %1/s
copars(54) =   0.1;         %1/s
copars(55) =       120;     %
copars(56) =    6;          %    
%Ke = 0;      %
%% NMDA
copars(57) =      5;        % 1/uM.s
copars(58) =      46.5;     % /s
copars(59) =    8.4;        % 1/s
copars(60) =    6.8;        % 1/s
copars(61) =     46.5;      % 1/s

copars(62) =    73.8;        % 1/s
copars(63) =    0.05;       % s
copars(64) =     0.005;     % s
copars(65) =     0.025;     % s
copars(66) =     0.003;     % s
copars(67) =     0.002;     %s
copars(68) =     10;        % mV
copars(69) =     -65;       % mV
copars(70) =      1;      % 
copars(71) =    40;     %mV

copars(72)=    -65.6e6;      %
copars(73)=    0.25*1e19/(2 * 1.602*6.022e23);%picoAmpere to Flux (micro mol/Ls)
copars(74)=     0.1;      %
copars(75)=          5;      %
copars(76)=             0.5;      %
copars(77)=            0.1;      %
copars(78)=   1 ;      %
copars(79)=    50;      %
copars(80)=   3 ;      %
copars(81)=     0.1;      %

%% AMPA
copars(82)= 15;      %
copars(83)= 1.8 ; % 1/uMs
copars(84)=  10;      %
copars(85)=  1.6e4;      %
copars(86)= 7e2;      %
copars(87)=  1e2;      %
copars(88)=  3e2;      %
copars(89)= 10;      %
copars(90)=  1.6e4;      %
copars(91)=  2.4e3;      %
copars(92)= 1e4;      %
copars(93)= 5e3;      %
copars(94)= 1.5e2;      %
copars(95)= 2.1;      %
copars(96)= 15;      %
copars(97)= 1e3;      %
copars(98)= 1.2e4;      %
%% Calmodulin

copars(99)= 0.001;      %
copars(100)= 01.0;      %
copars(101)=0.01;      %
copars(102)=0.01;      %
copars(103)=0.005;      %
copars(104)= 1;      %


end

