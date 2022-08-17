figure(412)
hold on
% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 1500;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% n = 2^nextpow2(L);
% Y = fft(y(:,pAMPK),n);
% Y(1) = [];
% P2 = abs(Y/L);
% P1 = P2(1:n/2+1);
% P1(:,2:end-1) = 2*P1(:,2:end-1);
% f = Fs*(0:(L/2))/L;
% plot(0:(Fs/n):(Fs/2-Fs/n),P1(1:n/2))
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

Y = fft(y(1:5000,12));
n = length(Y);
power = abs(Y(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 100;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
plot(freq,power)
xlabel('Freq')
ylabel('Power')