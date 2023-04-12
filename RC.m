function [outpRC, outpRRC] = RC(inp,Tsymbol,N,Fs)
%RC Summary of this function goes here
%   Detailed explanation goes here
%   N amount of taps for the filter

b = 0.3;
fc = 1e6;  % cutoff freq

% N = length(inp);  % moet groot genoeg zijn zodat je genoeg van sinc hebt
% Fs = 3e6;   % sampling freq > 2*fc
Ts = 1/Fs;
t = 0:1/Fs:(N-1)*Ts ;
delta_f = Fs/N;
fmax = delta_f*(N)/2;
fvector = linspace(-fmax,fmax,N+1);
fvector = fvector(1:end-1);
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

% figure
% plot(H)
% title('spectrum of filter')


h_rc = ifftshift(ifft(fftshift(H)));
ampl = max(h_rc);

h_rc = h_rc/ampl;
H_rrc = sqrt(H/ampl);
h_rrc = ifftshift(ifft((H_rrc)));

% figure
% plot(t,h_rc)
% x = 0:Tsymbol:(N-1)*Ts;
% hold on
% plot(x,zeros(1,length(x)),'x')
% 
% title('impulse response RC')

% figure
% plot(t,h_rrc)
% x = 0:Tsymbol:(N-1)*Ts;
% hold on
% plot(x,zeros(1,length(x)),'x')
% 
% title('impulse response RRC')

outpRC = h_rc;
outpRRC = h_rrc;
% % figure
% % plot(H)
% % title('spectrum of filter')
% 
%  % spectrum shifted
% Hshift = ifftshift(H);
% Hrcshift = ifftshift(Hrc);
% % figure
% % plot(Hshift)
% % title('spectrum of filter shifted')
% 
% 
% % ifft van dat spectrum
% ifftHshift = ifft(Hshift);
% ifftHrcshift = ifft(Hrcshift);
% % figure
% % plot(t,ifftHshift)
% % title('time signal of shifted spectrum')
% 
% 
% % terug shiften in tijdsdomein + normalisatie
% Hend = fftshift(ifftHshift/ifftHshift(1));
% Hrcend = fftshift(ifftHrcshift/ifftHrcshift(1));
% % Hrcend = sqrt(Hend);
% figure
% plot(t,Hend)
% 
% x = 0:Tsymbol:(N-1)*Ts;
% 
% hold on
% plot(x,zeros(1,length(x)),'x')
% 
% title('impulse response RC')
% 
% figure
% plot(t,Hrcend)
% hold on
% plot(x,zeros(1,length(x)),'x')
% 
% title('impulse response RRC')
% % outpRC = Hend;
% % outpRRC = Hrcend;
% % close all

end

% filtered = [];
% for i=1:length(inpf)
% 
%     f = (Fs/N * (i-1)) - Fs/2; 
% 
%     if  abs(f) <= (1-b)/(2*Tsymbol)
% 
%         filtered = [filtered inpf(i).*Tsymbol];
% 
%     elseif (1-b)/(2*Tsymbol) <= abs(f) && abs(f) <= (1+b)/(2*Tsymbol)
% 
%         filtered = [filtered inpf(i) * Tsymbol/2*(1 + cos( (pi*Tsymbol/b)*(abs(f) - ((1-b)/(2*Tsymbol)) )))];
% 
%     elseif (1+b)/(2*Tsymbol) < abs(f)
% 
%         filtered =  [filtered 0];
% 
%     end

% Real check
% ISI
% 
% figure
% plot(filtered)
% title('spectrum of filter')
% 
%  % spectrum shifted
% filtershifted = ifftshift(filtered);
% figure
% plot(filtershifted)
% title('spectrum of filter shifted')
% 
% 
% % ifft van dat spectrum
% ifftfilter = ifft(filtershifted);
% figure
% plot(t,ifftfilter)
% title('time signal of shifted spectrum')
% 
% 
% % terug shiften in tijdsdomein + normalisatie
% ifiltershifted = fftshift(ifftfilter/ifftfilter(1));
% figure
% plot(t,ifiltershifted)
% 
% x = 0:Tsymbol:(N-1)*Ts;
% 
% 
% hold on
% plot(x,zeros(1,length(x)),'x')
% 
% title('impulse response after shift')
%  
% outp = ifiltershifted;