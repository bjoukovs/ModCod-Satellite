%% Same main as test_mapping_demapping_filter_v2 but with variation of the number of iterations of hard decoding

clear all;
close all;
addpath('LDPC');

format long

fsymbol=1e6;
T=1/fsymbol;

U=4;
fs = fsymbol*U;
ts = 1/fs;

M = 16;
bits_per_symbol = log2(M)

blocklength=128;
bits = randi(2,bits_per_symbol*100*blocklength,1); %100k symbols
% ATTENTION put more than bits_per_symbol*100*blocklength
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;

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
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(encoded_message,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=4;
%figure;
%stem(symb_tx)
symb_tx = upsample(symb_tx,U);
%figure;
%stem(symb_tx)

%BER_moyen=[];

%% Computing different parts of H with removed columns
[row,col]=size(H);
subH_cell=cell(row,col);
for i=1:row
    for j=1:col
        if H(i,j)==1
            %calculate syndrome WITHOUT the validation node j (tip:
            %we do it only for the check node i to save time.)
            %Don't forget to discard the validation node j !
            subH_cell{i,j} = horzcat(H(i,1:j-1),H(i,j+1:end));
        end
    end
end

nbre_iter_hard=[1 3 5 10 20 30];
EbN0 = logspace(0,2.5,25);
BER = zeros(length(EbN0),length(nbre_iter_hard)); %the pth column contains the BER for the pth number of iteration of hard decoding
%% Loop with different hard decoding nbre of iterations
for p=1:length(nbre_iter_hard)
    p
    nbre_iter_hard(p)
    % Loop for different bit energies +  calculating BER
    %EbN0 = logspace(0,2,5);
    %EbN0 = [10^1.5]

    %EbN0 = logspace(0,8,5);
    %EbN0 = linspace(0,100,100)
    %EbN0=1:1:100;
    %BER = zeros(length(EbN0),1);
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
        bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);

        %%%%%%%%%%%%%%% CHANNEL DECODING %%%%%%%%%%%%%%%%%%%%
        % Now we have to cut in blocks of 2*128 bits to take the redundancy
        bits_rx_dec=[];

        for j=1:length(bits_rx)/(2*blocklength)
            j
            received = hard_decoderv3(H,bits_rx((j-1)*blocklength*2+1:j*blocklength*2,1)',nbre_iter_hard(1,p),subH_cell);
            received = received';
            bits_rx_dec=[bits_rx_dec;received(blocklength+1:end)]; %to take only the information bits (not the redudancy)
            %here we take only the second part of each decoded block since it
            %contains the message without redundancy
        end



        % Check error
        BER(i,p) = bit_error_rate(bits_rx_dec, bits);
        %BER(i,2) = bit_error_rate(encoded_message, bits_rx);
        %BER_moyen(i)=BER_moyen(i)+bit_error_rate(bits, bits_rx);
        %end
        %BER_moyen(i)=BER_moyen(i)/5;
        i
    end

    figure;
    semilogy(10*log10(EbN0),BER(:,p));hold on;
    %figure(25);semilogy(10*log10(EbN0),BER_moyen);
    xlabel("Eb/N0 (dB)"); hold on;
    ylabel("Bit error rate"); hold on;
end