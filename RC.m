function [outp] = RC(inp,Tsymbol)
%RRC Summary of this function goes here
%   Detailed explanation goes here

b = 0.3;
fc = 1e6;  % cutoff freq

N = length(inp);  % moet groot genoeg zijn zodat je genoeg van sinc hebt
Fs = 3e6;   % sampling freq > 2*fc
Ts = 1/Fs;
t = 0:1/Fs:(N-1)*Ts ;
inpf = abs(fft(inp,N));  %fft

filtered = [];

for i=1:length(inpf)

    f = (Fs/N * (i-1)) - Fs/2; 

    if  abs(f) <= (1-b)/(2*Tsymbol)

        filtered = [filtered inpf(i).*Tsymbol];

    elseif (1-b)/(2*Tsymbol) <= abs(f) && abs(f) <= (1+b)/(2*Tsymbol)

        filtered = [filtered inpf(i) .* Tsymbol/2*((1+cos((pi*Tsymbol/b)*(abs(f)-(1-b)/(2*Tsymbol)))))];

    elseif (1+b)/(2*Tsymbol) < abs(f)

        filtered =  [filtered 0];

    end
end

% Real check
% ISI



figure
plot(filtered)
title('spectrum of filter')

 % spectrum shifted
filtershifted = fftshift(filtered);
figure
plot(filtershifted)
title('spectrum of filter shifted')


% ifft van dat spectrum
ifftfilter = ifft(filtershifted);
figure
plot(t,ifftfilter)
title('time signal of shifted spectrum')


% terug shiften in tijdsdomein
ifiltershifted = ifftshift(ifftfilter);
figure
plot(t,ifiltershifted)

x = 0:Ts:(N-1)*Ts;

hold on
plot(x,zeros(1,length(x)),'x')

title('impulse response after shift')




outp = ifiltershifted;

end
