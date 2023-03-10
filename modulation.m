clc
clear
close all

%% PAM

bits= [1 0 0 1 1 0 1 0 1 0 0 1]';

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

bits= [1 0 0 1 1 0 1 0 1 0 0 1 1 1 0 1]';

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

demappedbits = demapping(mappedbits,Nbps,'qam');

figure
stem(demappedbits)
title(['mapped bitstream after QAM-' num2str(2^Nbps) ' demapping'])

error = norm(bits - demappedbits);



%% Rectangular window
clc
clear
close all

%tx
bits= [1 0 1 1 0 0 1 0]'
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

Tsymbol = 1e-6;  % symbol length
M = 5;

%bit stream
bits= randi([0 1], 1,8)'

%QAM-16 modulation
mappedbits = mapping(bits,4,'qam')

% Upsample function
tx =  sampler(mappedbits,M,'up');

% apply RRC 
ty = RRC(tx,Tsymbol)

%%% Channel %%%

% Apply AWGN
%y = awgn(x,20,'measured');  % SNR = 20dB idk wat measured voor staat
% zelf implementeren

% Rx %

% Apply RRC again
ry = RRC(ty,Tsymbol);

% Downsampling
rx =  sampler(ry,M,'down')';

% QAM-16 demodulation
receivedbits = demapping(rx,4,'qam')

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

Tsymbol = 1/(1.5e6);

test = zeros(1,101);
test(51) = 1;

RRC(test,Tsymbol)


% Vragen: wat is Tsymbol en Tsample nu en heo verhouden die zich to elkaar






