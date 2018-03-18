clear all;
close all;

M = 16;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*100,1);
bits = bits -1;

%% MAPPING
modulation = "pam";
%modulation = "qam";

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% Filter
cutf = 1e6; %cutoff frequency
fs = 2e6;
symb_tx = halfrootfilter(symb_tx,cutf,0,fs);

%% DOWNSAMPLING
M=2;
symb_tx = downsample(symb_tx,M);
figure;stem(symb_tx);



%% DEMAPPING

bits_rx = demapping(symb_tx,bits_per_symbol,modulation)

%% Check error

error_norm = norm(bits_rx - bits)