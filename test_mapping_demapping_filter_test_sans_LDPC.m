clear all;
close all;
addpath('LDPC');

format long

fsymbol=1e6;
T=1/fsymbol;

U=4;
fs = fsymbol*U;
ts = 1/fs;

M = 64;
bits_per_symbol = log2(M)

%blocklength=32;
bits = randi(2,bits_per_symbol*100000,1); %100k symbols
% ATTENTION put more than bits_per_symbol*100*blocklength
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;

bits = bits -1;


%% MAPPING
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=4;
%figure;
%stem(symb_tx)
symb_tx = upsample(symb_tx,U);
%figure;
%stem(symb_tx)

%BER_moyen=[];

%% Loop for different bit energies +  calculating BER
EbN0 = logspace(-0.4,2,10);

%EbN0 = logspace(0,8,5);
%EbN0 = linspace(0,100,100)
%EbN0=1:1:100;
BER = zeros(length(EbN0),1);
for i=1:length(EbN0)
    %BER_moyen(i)=0;
    %for j=1:10
    % First half root filter
    beta = 0.3; %imposed
    symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs);
    %figure;
    %stem(symb_tx_filtered)
    
    
    % Adding noise
    symb_tx_noisy = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));
    % length(encoded_message) or length(bits) in last argument of AWNG?
    
    %figure;
    %stem(symb_tx_noisy)
    
    
    % Second half root filter
    beta = 0.3; %imposed
    
    symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs);
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
    %BER_moyen(i)=BER_moyen(i)+bit_error_rate(bits, bits_rx);
    %end
    %BER_moyen(i)=BER_moyen(i)/5;
    i
end

figure;
semilogy(10*log10(EbN0),BER);
%figure(25);semilogy(10*log10(EbN0),BER_moyen);
xlabel("Eb/N0 (dB)");
ylabel("Bit error rate");