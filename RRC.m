function [outpRRC] = RRC(Tsymbol,N,Fs,b)
%RC Summary of this function goes here
%   Detailed explanation goes here
%   N amount of taps for the filter
%   Fs = sampling freq > 2*fc 

delta_f = Fs/N;
fmax = delta_f*(N-1)/2;
fvector = linspace(-fmax,fmax,N);
H = zeros(1,N);

for i=1:N
    f = fvector(i);

    if  abs(f) < (1-b)/(2*Tsymbol)
        H(i) = Tsymbol;

    elseif (1-b)/(2*Tsymbol) <= abs(f) && abs(f) <= (1+b)/(2*Tsymbol)
        H(i) = Tsymbol/2*(1 + cos( (pi*Tsymbol/b)*(abs(f) - ((1-b)/(2*Tsymbol)) )));

    elseif (1+b)/(2*Tsymbol) < abs(f)
        H(i) = 0;
    end
end




Hrrc = sqrt(H);
outpRRC = fftshift(ifft(ifftshift(Hrrc)));

%Normalize
amplitude = sqrt(max(real(conv(outpRRC, outpRRC))));
outpRRC = outpRRC./amplitude;

%Plots
Ts = 1/Fs;
t = 0:1/Fs:(N-1)*Ts ;

figure
plot(H)
title('Frequency response Raised Cosine filter')

figure
H_RRC = fftshift(ifft(ifftshift(H./amplitude^2)));
plot(t,H_RRC,'LineWidth',1);
x = 0:Tsymbol:(N-1)*Ts;
middle = x((end+1)/2);
x((end+1)/2) = [];
hold on
plot(t,outpRRC,'LineWidth',1);
hold on
plot(x,zeros(1,length(x)),'x','LineWidth',2)
hold on
plot(middle,1,'x','LineWidth',2);
% plot(t,H_RRC,'x','MarkerIndices',1:1/Fs:length(x)+1)
title('Impulse response Raised Cosine filter')
legend('Raised Cosine','Root Raised Cosine','Sampled every Tsymbol')


end