clc
clear
close all

%% PAM

bits= [1 0 0 1 1 0 1 0 1 0 0 1].';

Nbps = 2;

figure
stem(bits)
title('bitstream before mapping')

mappedbits = mapping(bits,Nbps,'pam');

figure
scatter(mappedbits,zeros(length(mappedbits)))
title('constellation diagram PAM')

figure
stem(mappedbits)
title(['bitstream after PAM-' num2str(Nbps + 2) ' mapping'])

demappedbits = demapping(mappedbits,Nbps,'pam');

figure
stem(demappedbits)
title(['mapped bitstream after PAM-' num2str(Nbps + 2) ' demapping'])

error = norm(bits - demappedbits);

%% QAM
clc
clear
close all

bits= [1 0 0 1 1 0 1 0 ].'

Nbps = 2;

figure
stem(bits)
title('bitstream before mapping')

mappedbits = mapping(bits,Nbps,'qam');

figure
plot(mappedbits,'x')
title('constellation diagram QAM')

figure
stem(mappedbits)
title(['bitstream after QAM-' num2str(2^Nbps) ' mapping'])

% % sampling
% x =  sampler(mappedbits,5,'up');
% y =  sampler(x,5,'down');

demappedbits = demapping(mappedbits,Nbps,'qam')
%demappedbits = demapping(y,Nbps,'qam')
figure
stem(demappedbits)
title(['mapped bitstream after QAM-' num2str(2^Nbps) ' demapping'])

error = norm(bits - demappedbits);



%% Rectangular window
clc
clear
close all

%tx
bits = [1 0 0 1 0 1 1 0].'
%QAM mod
mappedbits = mapping(bits,4,'qam');
%Upsample
x =  sampler(mappedbits,10,'up');
%Filter
window = rectwin(10);
idk = conv(mappedbits,window);

%rx
%Filter
idkrx = conv(idk,window);
%Downsample
x =  sampler(idkrx,10,'down');
%QAM demod
receivedbits = demapping(x,4,'qam')


%% Root raised cosine
% Tx %
clc
clear
close all


Tsymbol = 1/(1e6);  % symbol length
M = 10;
Nbps = 2;           %number of bits per symbol
Fs = 3e6;           %sampling frequency
% EbN0 = 5;           %energy per bit to noise power spectral density ratio
EbN0 = 1:0.05:20;
finalbits = [];

bits = randi([0 1], 1,20).'

for i=1:length(EbN0)


%bit stream

%bits = [1 0 0 1 0 1 0 0].'

%QAM mod
mappedbits = mapping(bits,Nbps,'qam');

%Upsample
x =  sampler(mappedbits,M,'up').';

%Filter

[~,filter] = RC(x,Tsymbol);

RRCtest = conv(x,filter,'same');

%Noise

RRCtestNoisy = AWGNoise(RRCtest,1/Fs,Nbps,EbN0(i));

%Filter

RRCy = conv(RRCtestNoisy,filter,'same');

%Downsample
y =  sampler(RRCy,M,'down');

%QAM demod
receivedbits = demapping(y,Nbps,'qam');

finalbits = [finalbits,receivedbits];
close all
end

numberOfErrors = [];
for i=1:length(EbN0)
    numberOfErrors = [numberOfErrors biterr(finalbits(:,i) ,bits)/length(bits)];
end


figure
plot(EbN0,numberOfErrors)
% figure
% stem(bits)
% hold on
% stem(receivedbits)
% 
% figure
% plot(mappedbits,'x')
% hold on
% plot(y,'.')
% figure
% plot(RRCy,'.')



%% sampler test

clear
clc

bits= [1 0 1 1 0 0 ];
x =  sampler(bits,3,'up');
y =  sampler(x,3,'down');

%% RC test

clear
clc
close all

Tsymbol = 1/(1e6);

test = zeros(20,1);
test(12) = 1;

[x,y] = RC(test,Tsymbol);

RCtest = conv(x,test);

RRCtest = conv(test,y);
RRCtest = conv(RRCtest,y,"same")/2.5;

figure
plot(RCtest)
hold on
plot(RRCtest)
%plot(circshift(RRCtest/2.54,-50))
title('lets see if it works')


