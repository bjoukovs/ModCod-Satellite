clear all;
close all;

fsymbol = 1e6;
T = 1/fsymbol;

fs = 2e6;
ts = 1/fs;


M = 16;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*1000,1);
bits = bits -1;

%% MAPPING
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=2;
symb_tx = repelem(symb_tx,U);

%% Rectangular Filter
% cutf = 1e6; %cutoff frequency
% fs = 2e6;
% symb_tx = halfrootfilter(symb_tx,cutf,0,fs,M);

%% Raised cosine filter
beta = 0.3; %imposed
fs = 2e6;
ts = 1/fs;

%symb_tx = halfroot_opti(symb_tx,U,beta,ts);

EbN0 = logspace(1,15);
BER = zeros(length(EbN0),1);
for i=1:length(EbN0)
    
    %% Adding noise
    symb_tx_noisy = AWNG(symb_tx,EbN0(1,i),M,U*fs,modulation);

    %% DOWNSAMPLING
    symb_tx_noisy = downsample(symb_tx_noisy,U);

    %% DEMAPPING

    bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);

    %% Check error

    BER(i) = bit_error_rate(bits, bits_rx);
end

loglog(EbN0,BER);