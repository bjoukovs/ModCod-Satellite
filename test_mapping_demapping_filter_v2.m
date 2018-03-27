clear all;
close all;

fsymbol=1e6;
T=1/fsymbol;

fs = 1e6;
ts = 1/fs;

M = 4;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*1000,1); %1k symbols
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;

SignalEnergy = (trapz(abs(bits).^2))*1/fs;
Eb = SignalEnergy/length(bits)/2;

bits = bits -1;

%% MAPPING
%modulation = 'pam';
modulation = "qam";

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=10;
symb_tx = upsample(symb_tx,U);

%% Loop for different bit energies +  calculating BER
EbN0 = logspace(1,8,100);
BER = zeros(length(EbN0),1);
for i=1:length(EbN0)
    % First half root filter
    beta = 0.3; %imposed
    symb_tx = halfroot_opti_v2(symb_tx,beta,T,U*fs); %ATTENTION fs ou U*Fs ?
    
    
    % Adding noise
    symb_tx_noisy = AWNG(symb_tx,EbN0(1,i),M,U*fs,modulation);
    
    
    % Second half root filter
    beta = 0.3; %imposed
    
    symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,U*fs);
    

    % DOWNSAMPLING
    symb_tx_noisy = downsample(symb_tx_noisy,U);
    

    % DEMAPPING
    bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);
    

    % Check error
    BER(i) = bit_error_rate(bits, bits_rx);
    return 
end

loglog(EbN0,BER);
xlabel("Eb/N0");
ylabel("Bit error rate");