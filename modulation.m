%% PAM
clc
clear
close all

bits= [1 0 0 1 1 0 1 0 1 0 0 1].';

Nbps = 1;

figure
stem(bits)
title('bitstream before mapping')

mappedbits = mapping(bits,Nbps,'pam');

figure
scatter(mappedbits,zeros(length(mappedbits)))
title('constellation diagram PAM')

figure
stem(mappedbits)
title(['bitstream after PAM-' num2str(Nbps *2) ' mapping'])

demappedbits = demapping(mappedbits,Nbps,'pam');

figure
stem(demappedbits)
title(['mapped bitstream after PAM-' num2str(Nbps*2) ' demapping'])

error = norm(bits - demappedbits,1)
BER = error/length(bits);

%% QAM
clc
clear
close all

bits = randi([0 1],1, 12000).';

Nbps = 6;

% figure
% stem(bits)
% title('bitstream before mapping')

mappedbits = mapping(bits,Nbps,'qam');

% figure
% plot(mappedbits,'o')
% title('constellation diagram QAM')
% 
% figure
% stem(mappedbits)
% title(['bitstream after QAM-' num2str(2^Nbps) ' mapping'])

demappedbits = demapping(mappedbits,Nbps,'qam');
% figure
% stem(demappedbits)
% title(['mapped bitstream after QAM-' num2str(2^Nbps) ' demapping'])

error = norm(bits - demappedbits,1);
BER = error/length(bits);

%% RRC test
clc
clear
close all

Fsymbol = 1e6;  % symbols/s
M = 10;
Fs = M*Fsymbol;   % sample freq
N = 1003;   % oneven -> real , even -> complex god weet wrm
Nbps = [1,2,4,6];
SNREb = 0:2:10;
Modulation = 'pam';
BERmatrix = zeros(length(Nbps),length(SNREb));

filter = RRC(1/Fsymbol,N,Fs,0.3);

% figure
%t = 0:1/Fs:(N-1)/Fs;
% plot(t,conv(filter,filter,"same"))
% hold on
% x = 0:1/Fsymbol:(N-1)/Fs;
% plot(x,zeros(length(x)),'*')
% title("RRC * RRC filter")

bits = randi([0 1],1, 1200000).';

tic
for i=1:length(Nbps)
     if Nbps(i) ~= 1
         Modulation = 'qam';
     end
    for j=1:length(SNREb)    
        mappedbits = mapping(bits,Nbps(i),Modulation);
        upmapped =  sampler(mappedbits,M,'up').';
        
        RRCtx = conv(upmapped,filter,"same");
        
        RRCtxnoisy = AWGNoise(RRCtx,Fs,Fsymbol,Nbps(i),SNREb(j),Modulation);
        
        RRCrx = conv(RRCtxnoisy,filter,"same");
        
        downmapped =  sampler(RRCrx,M,'down');
%          figure
%          plot(mappedbits,'o')
%          hold on
%          plot(downmapped,'*')
%          title('constellation diagram QAM')
%          legend('Transmitted symbols','Received symbols')

        demapped = demapping(downmapped,Nbps(i),Modulation); % 2*600
%         demapped = reshape(demapped,size(bits)); % 1*1200
        
        error = norm(bits - demapped,1);
        BER = error/length(bits);
        BERmatrix(i,j) = BER;
    end
end

figure
for k=1:height(BERmatrix)
    plot(SNREb,BERmatrix(k,:),'-o')
    hold on
end

toc

hold off
title('BER in function of Eb/N0')
ylabel('BER')
xlabel('Eb/N0 [dB]')
legend('BPSK','QPSK','QAM-16','QAM-64')
set(gca, 'YScale', 'log')


%% LDPC step 1
clc
clear
close all

H = [1 1 0 1 1 0 0 1 0 0;
     0 1 1 0 1 1 1 0 0 0;
     0 0 0 1 0 0 0 1 1 1;
     1 1 0 0 0 1 1 0 1 0;
     0 0 1 0 0 1 0 1 0 1];

 
H1 = mod(rref(H), 2);
K = 5;
M = 10;
G = [H1(:, M - K + 1 : M)', eye(K)];
bits = [1 1 0 1 0];
enc_bits_og = LDPC_enc(G,bits);

%Random error
error = zeros(1,length(enc_bits_og));
rdidx = randi([1 length(enc_bits_og)],1);
error(rdidx) = 1;

enc_bits = xor(enc_bits_og,error);
enc_bits = double(enc_bits);

% enc_bits = enc_bits_og; %uncomment to remove the error


dec_bits = LDPC_harddecode(enc_bits,H1,100);
enc_bits_test = enc_bits_og(6:end);

if(isequal(enc_bits_test, dec_bits))
    disp('Decoded succesfully!')
else
    disp('Decoding failed!')
end



%% LDPC encoder part 2 (hard decoding)
clc
clear
close all
% ========================================================================
% Graph for every modulation type, for every modulation evaluate BER /
% Eb/N0 for several number of max iterations
Fsymbol = 1e6;  % symbols/s
M = 10;
Fs = M*Fsymbol;   % sample freq
N = 1003;   % oneven -> real , even -> complex god weet wrm
Nbps = [1,2,4,6];
SNREb = -2:1:10;
Modulation = 'pam';
BERmatrix = zeros(length(Nbps),length(SNREb)); %----------

filter = RRC(1/Fsymbol,N,Fs,0.3);

% Encode bits
BlockSize = 128;
n_Blocks = 128;
code_rate = 1/2;
og_bits = randi([0 1],1, BlockSize*n_Blocks).';
bits = og_bits;
H0 = generate_ldpc(BlockSize, BlockSize / code_rate,0,1,3); % Create initial parity check matrix of size 128 x 256
CodeBlocks = reshape(bits,BlockSize,n_Blocks); % Every row a block of 128 bits
[parity_bits, H] = encode_ldpc(CodeBlocks, H0, 0); % Compute parity bits (128 x number of packets) and generate final parity check matrix
encoded_block = [parity_bits; CodeBlocks];
encoded_bits = reshape(encoded_block,BlockSize*n_Blocks/ code_rate,1);

n_iterations = [0 1 2 5 10 20];
BERmatrix = zeros(length(n_iterations),length(SNREb));

for i=1:length(n_iterations)
    if n_iterations(i) ~= 0
        bits = encoded_bits;
    end
    for j=1:length(SNREb)    
        mappedbits = mapping(bits,Nbps(1),Modulation);
        upmapped = sampler(mappedbits,M,'up').';

        RRCtx = conv(upmapped,filter,"same");

        RRCtxnoisy = AWGNoise(RRCtx,Fs,Fsymbol,Nbps(1),SNREb(j),Modulation);

        RRCrx = conv(RRCtxnoisy,filter,"same");

        downmapped =  sampler(RRCrx,M,'down');

        demapped = demapping(downmapped,Nbps(1),Modulation);
        
        if n_iterations(i) ~= 0
           demapped = reshape(demapped, BlockSize/code_rate,n_Blocks);
           
           % Decoding each block
           dec_bits = zeros(BlockSize,n_Blocks);
           for k = 1:n_Blocks
               dec_block = LDPC_harddecode(demapped(:,k)',H,n_iterations(i));
               dec_bits(:,k) = dec_block';
           end
           
           dec_bits = reshape(dec_bits,1,BlockSize*n_Blocks);
        else 
           dec_bits = demapped;
        end
        
        if n_iterations(i) ~= 0
            error = norm(og_bits - dec_bits',1);
        else
            error = norm(og_bits - dec_bits,1);
        end
        
        BER = error/length(og_bits);
        BERmatrix(i,j) = BER;
    end
end

figure
for k=1:height(BERmatrix)
    plot(SNREb,BERmatrix(k,:),'-')
    hold on
end

hold off
title('BER in function of Eb/N0 + Hard decoding BPSK')
ylabel('BER')
xlabel('Eb/N0 [dB]')
legend('No encoding','1 iteration','2 iterations','5 iterations','10 iterations','20 iterations')
set(gca, 'YScale', 'log')


%% LDPC soft decoder (pt3)















