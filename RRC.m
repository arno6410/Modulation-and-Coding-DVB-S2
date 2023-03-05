function [outp] = RRC(inp,Tsymbol)
%RRC Summary of this function goes here
%   Detailed explanation goes here

b = 0.3;
fc = 1e6;  % cutoff freq

N = length(inp);  % moet groot genoeg zijn zodat je genoeg van sinc hebt
Fs = 3e6;   % sampling freq > 2*fc
inpf = abs(fft(inp,N));  %fft

filtered = [];

for i=1:length(inpf)

    f = (Fs/N * (i-1)) - Fs/2; 

    if  abs(f) <= (1-b)/(2*Tsymbol)

        filtered = [filtered sqrt(inpf(i).*Tsymbol)];

    elseif (1-b)/(2*Tsymbol) <= abs(f) && abs(f) <= (1+b)/(2*Tsymbol)

        filtered = [filtered inpf(i) .* sqrt(Tsymbol/2*(1/2 *(1+cos((pi*Tsymbol/b)*(abs(f)-(1-b)/(2*Tsymbol))))))];

    elseif (1+b)/(2*Tsymbol) < abs(f)

        filtered =  [filtered 0];

    end
end

outp = ifft(filtered);

end

