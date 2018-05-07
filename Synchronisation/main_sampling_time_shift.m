%% Same main as test_mapping_demapping_filter_v2 but with sampling time shift included
%% ATTENTION here no channel coding
%% ATTENTION DOES NOT WORK FOR PAM FOR THE MOMENT
%% because PAM decoder does not handle complex numbers

clear all;
close all;
addpath('..'); %parent directory

format long

fsymbol=1e6;
T=1/fsymbol;

U=4;
fs = fsymbol*U;
ts = 1/fs;

M = 16;
bits_per_symbol = log2(M)

blocklength=128;
bits = randi(2,bits_per_symbol*10*blocklength,1); %100k symbols
% ATTENTION put more than bits_per_symbol*100*blocklength
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;

bits = bits -1;


%%%% CHANNEL CODING %%%%%
% encoded_message=[];
% H0 = makeLdpc(blocklength,2*blocklength,0,1,3);   %128*256 H -> encodes 128 bits
% for i=1:length(bits)/blocklength %we divide the bits in blocks of 128
%     i
%     [codedbits, H] = makeParityChk(bits((i-1)*blocklength+1:i*blocklength,1), H0, 0); 
%     encoded_message = [encoded_message;codedbits;bits((i-1)*blocklength+1:i*blocklength,1)];
% end


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

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
U=100;
symb_tx_upsampled = upsample(symb_tx,U);

t_original = linspace(0,length(symb_tx_upsampled)*ts/U,length(symb_tx_upsampled));
n_original = 1:U:length(symb_tx_upsampled);

sample_time_shift = linspace(0,0.05,10); %expressed in the form of a percentage of the original sampling time


%% Computing different parts of H with removed columns
% [row,col]=size(H);
% subH_cell=cell(row,col);
% for i=1:row
%     for j=1:col
%         if H(i,j)==1
%             %calculate syndrome WITHOUT the validation node j (tip:
%             %we do it only for the check node i to save time.)
%             %Don't forget to discard the validation node j !
%             subH_cell{i,j} = horzcat(H(i,1:j-1),H(i,j+1:end));
%         end
%     end
% end



%% Loop for different bit energies +  calculating BER
EbN0 = logspace(0,2,10);
%EbN0 = [10^1.5]

%EbN0 = logspace(0,8,5);
%EbN0 = linspace(0,100,100)
%EbN0=1:1:100;

BER = zeros(length(EbN0),length(sample_time_shift));

t=[0:ts:(length(symb_tx)-1)*ts]';


for p=1:length(sample_time_shift)
    
    T_sampling = ts/U*bits_per_symbol;
    n_sampling = n_original + ones(1,length(n_original))*U*sample_time_shift(p);
    
    for i=1:length(EbN0)
        
        % First half root filter
        beta = 0.3; %imposed
        symb_tx_filtered = halfroot_opti_v2(symb_tx_upsampled,beta,T,fs);

        % Adding noise
        [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));

        % Second half root filter
        beta = 0.3; %imposed

        symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs);



        % DOWNSAMPLING TAKING INTO ACCOUNT THE SCO
        sampled_signal = zeros(length(symb_tx),1);
        
        for j=1:length(n_sampling)
           if ceil(n_sampling(j)) <= length(symb_tx_noisy)
               
               n_low = floor(n_sampling(j));
               n_high = ceil(n_sampling(j));
               
               %case if we fall just right on a sample
               if n_low == n_high
                   sampled_signal(j) = symb_tx_noisy(n_low);
               else
                   %we have to interpolate linearly
                   s_low = symb_tx_noisy(n_low);
                   s_high = symb_tx_noisy(n_high);
                   
                   sampled_signal(j) = s_low + (s_high-s_low) * (n_sampling(j)-n_low);
               end
           end
        end
        

        % DEMAPPING
        bits_rx = demapping(sampled_signal,bits_per_symbol,modulation);

        %%%%%%%%%%%%%%% CHANNEL DECODING %%%%%%%%%%%%%%%%%%%%
        % Now we have to cut in blocks of 2*128 bits to take the redundancy
%         bits_rx_dec=[];
% 
%         for j=1:length(bits_rx)/(2*blocklength)
%             j
%             received = hard_decoderv3(H,bits_rx((j-1)*blocklength*2+1:j*blocklength*2,1)',10,subH_cell);
%             received = received';
%             bits_rx_dec=[bits_rx_dec;received(blocklength+1:end)]; %to take only the information bits (not the redudancy)
%             %here we take only the second part of each decoded block since it
%             %contains the message without redundancy
%         end



        % Check error
        BER(i,p) = bit_error_rate(bits_rx, bits);
        %BER(i,2) = bit_error_rate(encoded_message, bits_rx);
        %BER_moyen(i)=BER_moyen(i)+bit_error_rate(bits, bits_rx);
        %end
        %BER_moyen(i)=BER_moyen(i)/5;
        i
    end
end

figure;
for p=1:length(sample_time_shift)
    semilogy(10*log10(EbN0),BER(:,p));
    %figure(25);semilogy(10*log10(EbN0),BER_moyen);
    xlabel("Eb/N0 (dB)");
    ylabel("Bit error rate");
    hold on;
end