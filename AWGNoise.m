function y = AWGNoise(x,Tsamp,Nbits,ratio)
%ratio = Eb/N0

% Calculate values
Eb = (Tsamp/2*Nbits) * sum(x);
N0 = Eb/ratio;

% Calculate noise
gaussnoiseRE = normrnd(0,1,size(x));
gaussnoiseIM = normrnd(0,1,size(x));
noise = sqrt(N0/Tsamp)*(gaussnoiseRE+1j*gaussnoiseIM);

% Add noise
y = x + noise;

end