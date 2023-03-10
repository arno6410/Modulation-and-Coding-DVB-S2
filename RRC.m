function [outp] = RRC(inp,Tsymbol)
%RRC Summary of this function goes here
%   Detailed explanation goes here

b = 0.3;
fc = 1e6;  % cutoff freq

N = length(inp);  % moet groot genoeg zijn zodat je genoeg van sinc hebt
Fs = 3e6;   % sampling freq > 2*fc
Ts = 1/Fs;
t = 0:1/Fs:(N-1)*Ts ;
inpf = abs(fft(inp,N));  %fft

H_RC = [];

for i=1:length(inpf)

    f = (Fs/N * (i-1)) - Fs/2; 

    if  abs(f) <= (1-b)/(2*Tsymbol)

        H_RC = [H_RC inpf(i).*Tsymbol];

    elseif (1-b)/(2*Tsymbol) <= abs(f) && abs(f) <= (1+b)/(2*Tsymbol)

        H_RC = [H_RC inpf(i) * Tsymbol/2*(1 + cos( (pi*Tsymbol/b)*(abs(f) - ((1-b)/(2*Tsymbol)) )))];

    elseif (1+b)/(2*Tsymbol) < abs(f)

        H_RC =  [H_RC 0];

    end
end

% Real check
% ISI



figure
plot(H_RC)
title('spectrum of filter')

 % spectrum shifted
H_RC = ifftshift(H_RC);
H_RRC = sqrt(H_RC);
figure
plot(H_RC)
title('spectrum of filter shifted')


% ifft van dat spectrum
h_RC = ifft(H_RC);
h_RRC = ifft(H_RRC);
figure
plot(t,h_RC)
title('time signal of shifted spectrum')


% terug shiften in tijdsdomein
h_RRC = fftshift(h_RRC/sqrt(h_RC(1)));
h_RC = fftshift(h_RC/h_RC(1)); 
figure
plot(t,h_RC)

x = 0:Tsymbol:(N-1)*Ts;


hold on
plot(x,zeros(1,length(x)),'.')

title('impulse response after shift')
 
 

outp = h_RC;

end

