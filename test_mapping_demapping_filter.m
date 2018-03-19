clear all;
close all;

M = 16;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*100,1);
bits = bits -1;

%% MAPPING
modulation = 'pam';
%modulation = "qam";

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol M times
M=2;
symb_tx = repelem(symb_tx,M);

%% Rectangular Filter
% cutf = 1e6; %cutoff frequency
% fs = 2e6;
% symb_tx = halfrootfilter(symb_tx,cutf,0,fs,M);

%% Raised cosine filter
beta = 0.3; %imposed
fs = 2e6;
ts = 1/fs;
figure(5)
stem(symb_tx)
symb_tx = halfroot_opti(symb_tx,M,beta,ts);

%% Adding noise
symb_tx_noisy = AWNG(symb_tx);

%% DOWNSAMPLING
symb_tx_noisy = downsample(symb_tx_noisy,M);
figure(6);stem(symb_tx_noisy);

%% DEMAPPING

bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation)

%% Check error

error_norm = norm(bits_rx - bits)