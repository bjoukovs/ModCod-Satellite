clear all;
close all;

fsymbol=1e6;
T=1/fsymbol;

fs = 2e6;
ts = 1/fs;

M = 16;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*100,1);

SignalEnergy = (trapz(abs(bits).^2))*1/fs;
Eb = SignalEnergy/length(bits)/2;

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

% figure(5)
% stem(symb_tx)
symb_tx = halfroot_opti_v2(symb_tx,beta,T,fs);

% %% Adding noise
% symb_tx_noisy = AWNG(symb_tx);
% 
% %% Filtering at reception
% symb_rx = halfroot_opti(symb_tx_noisy,M,beta,ts);
% %% DOWNSAMPLING
% symb_rx = downsample(symb_rx,M);
% figure(6);stem(symb_rx);
% 
% %% DEMAPPING
% 
% bits_rx = demapping(symb_rx,bits_per_symbol,modulation);
% 
% %% Check error
% 
% error_norm = norm(bits_rx - bits);