function y = AWGNoise(x,Fs,Fsymbol,Nbps,ratio,modulation)
%ratio = Eb/N0 in dB!!!

power_signal = rms(x)^2;
Eb = power_signal/(Nbps*Fsymbol);  % Energy per bit
N0 = Eb/10^(ratio/10);

Re_noise = randn(length(x),1)*sqrt(N0*Fs);

switch modulation
    case 'pam'
    y = x + Re_noise;
    case 'qam'
    Im_noise = randn(length(x),1)*sqrt(N0*Fs);
    y = x + Re_noise + 1j*Im_noise;
end

end