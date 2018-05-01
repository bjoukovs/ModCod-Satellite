clear all;
close all;
addpath('LDPC');

format long

fsymbol=1e6;
T=1/fsymbol;

U=4;
fs = fsymbol*U;
ts = 1/fs;

M = 2;
bits_per_symbol = log2(M)

blocklength=32;
bits = randi(2,bits_per_symbol*4000*blocklength,1); %100k symbols

bits = bits -1;

%%%% CHANNEL CODING %%%%%
encoded_message=[];
H0 = makeLdpc(blocklength,2*blocklength,0,1,3);   %128*256 H -> encodes 128 bits
for i=1:length(bits)/blocklength %we divide the bits in blocks of 128
    i
    [codedbits, H] = makeParityChk(bits((i-1)*blocklength+1:i*blocklength,1), H0, 0); 
    encoded_message = [encoded_message;codedbits;bits((i-1)*blocklength+1:i*blocklength,1)];
end
% We divided the message in blocks of 128 bits, added 128 bits of
% redundancy, thus we end with a sequence encoded_message that is twice as
% long as the bits

% On peut tout calculer en une fois en mettant les séquences de bits en
% colonnes dans une matrice (reshape)

%ATTENTION H remains the same
%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAPPING
modulation = 'pam';
%modulation = 'qam';

symb_tx = mapping(encoded_message,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=4;
%figure;
%stem(symb_tx)
symb_tx = upsample(symb_tx,U);
%figure;
%stem(symb_tx)



%% Loop for different bit energies +  calculating BER
EbN0 = logspace(-0.4,1.4,18);

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
    [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(encoded_message));
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
    %bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);
    
    %%%%%%%%%%%%%%% CHANNEL DECODING %%%%%%%%%%%%%%%%%%%%
    % Now we have to cut in blocks of 2*128 bits to take the redundancy
    bits_rx_dec=[];

    for j=1:length(symb_tx_noisy)/(2*blocklength)
        j
        received = soft_decoder_log(H,symb_tx_noisy((j-1)*blocklength*2+1:j*blocklength*2,1)',35,sigma);
        
        bits_rx_dec=[bits_rx_dec;received(blocklength+1:end,1)]; %to take only the information bits (not the redudancy)
        %here we take only the second part of each decoded block since it
        %contains the message without redundancy
    end
    
    

    % Check error
    BER(i) = bit_error_rate(bits, bits_rx_dec);
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