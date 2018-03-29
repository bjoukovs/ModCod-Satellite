clear all;
close all;

format long

fsymbol=1e6;
T=1/fsymbol;

fs = 1e6;
ts = 1/fs;

M = 16;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*10000,1); %1k symbols
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;

SignalEnergy = (trapz(abs(bits).^2));
Eb = SignalEnergy/length(bits)/2;

bits = bits -1;

%% MAPPING
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=5;
%figure;
%stem(symb_tx)
symb_tx = upsample(symb_tx,U);
%figure;
%stem(symb_tx)

%% Loop for different bit energies +  calculating BER
EbN0 = logspace(-0.4,2,100);
%EbN0 = linspace(0,100,100)
%EbN0=1:1:100;
BER = zeros(length(EbN0),1);
for i=1:length(EbN0)
    % First half root filter
    beta = 0.3; %imposed
    symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,U*fs); %ATTENTION fs ou U*Fs ?
    %figure;
    %stem(symb_tx_filtered)
    
    
    % Adding noise
    symb_tx_noisy = AWNG(symb_tx_filtered,EbN0(1,i),M,U*fs,modulation);
    %figure;
    %stem(symb_tx_noisy)
    
    
    % Second half root filter
    beta = 0.3; %imposed
    
    symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,U*fs);
    %figure;
    %stem(symb_tx_noisy)
    

    % DOWNSAMPLING
    symb_tx_noisy = downsample(symb_tx_noisy,U);
    %figure;
    %stem(symb_tx_noisy)
    

    % DEMAPPING
    bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);
    

    % Check error
    BER(i) = bit_error_rate(bits, bits_rx);
     
end

figure;
semilogy(10*log10(EbN0),BER);
xlabel("Eb/N0 (dB)");
ylabel("Bit error rate");