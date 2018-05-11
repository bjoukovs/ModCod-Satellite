%% Same main as test_mapping_demapping_filter_v2 but with CFO included
%% ATTENTION here no channel coding

clear all;
close all;
addpath('LDPC');
addpath('..'); %parent directory

format long

fsymbol=1e6;
T=1/fsymbol;

U=100;
fs = fsymbol*U;
ts = 1/fs;

M = 16;
bits_per_symbol = log2(M)

SIZE_PILOT=4*bits_per_symbol; %nbre of symbols in pilot
LENGTH_FRAME=300-SIZE_PILOT;

%blocklength=128;
bits = randi(2,bits_per_symbol*10*LENGTH_FRAME,1); %100k symbols
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

% On peut tout calculer en une fois en mettant les s�quences de bits en
% colonnes dans une matrice (reshape)

%ATTENTION H remains the same
%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAPPING
modulation = 'pam';
%modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% INSERTING PILOTS

pilot=makePilot(modulation,bits_per_symbol,SIZE_PILOT);
symb_tx=intoFrames(symb_tx,LENGTH_FRAME,pilot); %attention now symb_tx is longer

bits_with_pilots=demapping(symb_tx,bits_per_symbol,modulation); %save the bits of the symbol containing pilots,
% if comparison is needed after

%% OVERSAMPLING = replicating each symbol U times
%figure;
%stem(symb_tx)
symb_tx = upsample(symb_tx,U);
%figure;
%stem(symb_tx)

%BER_moyen=[];

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
%EbN0 = logspace(0,2,10);
EbN0 = 100;
%EbN0 = [10^1.5]

%EbN0 = logspace(0,8,5);
%EbN0 = linspace(0,100,100)
%EbN0=1:1:100;

%%%%% To investigate CFO only %%%%%
f_carrier=2e9; %2GHz
CFO=0:10^(-6)*f_carrier:10*10^(-6)*f_carrier;
phi0=0;
delta_f_tild=zeros(1,length(CFO));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% To investigate phi0 only %%%%%
% CFO=0;
% phi0=0:10:90;
% phi0=deg2rad(phi0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BER = zeros(length(EbN0),length(CFO));

%for p=1:length(phi0)
for p=1:length(CFO)
    p
    for i=1:length(EbN0)
        %BER_moyen(i)=0;
        %for j=1:10
        % First half root filter
        beta = 0.3; %imposed
        symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs,U);
        %figure;
        %stem(symb_tx_filtered)


        % Adding noise
        [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));


        %% Reception: multiplying with exp(j*CFO*t+phi0) to take the effect CFO and phase offset (cf. slide 10)
        %%%%% To investigate CFO only %%%%%
        
        RRCTaps=25*U+1;
        t=[-(RRCTaps/2)*ts : ts : ((length(symb_tx_noisy)-1)-RRCTaps/2)*ts]';
        symb_tx_noisy = symb_tx_noisy.*exp(1j.*(CFO(1,p).*t+phi0));


        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
         % Second half root filter
         beta = 0.3; %imposed 
         symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs,U);
         %%%% Removing the extra samples due to the 2 convolutions %%%%
         symb_tx_noisy=symb_tx_noisy(RRCTaps:end-RRCTaps+1);

%         % DOWNSAMPLING
          symb_tx_noisy = downsample(symb_tx_noisy,U);

          
        %%%% CFO estimation: we estimate the CFO on one frame %%%%
        %%%% There is no need to compensate it because it is too
        %%%% complicated
        check_length=SIZE_PILOT+LENGTH_FRAME;
        K=8;
        b1=1;
        b2=check_length+1;
        delta_f_tild(1,p) = diff_corr(symb_tx_noisy(b1:b2),pilot,K,CFO(1,p),T);
        
          
          
%         % DEMAPPING
%         if strcmp(modulation,'pam')
%             symb_tx_noisy=real(symb_tx_noisy);
%         end
%         bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);
% 
%         %%%%%%%%%%%%%%% CHANNEL DECODING %%%%%%%%%%%%%%%%%%%%
%         % Now we have to cut in blocks of 2*128 bits to take the redundancy
% %         bits_rx_dec=[];
% % 
% %         for j=1:length(bits_rx)/(2*blocklength)
% %             j
% %             received = hard_decoderv3(H,bits_rx((j-1)*blocklength*2+1:j*blocklength*2,1)',10,subH_cell);
% %             received = received';
% %             bits_rx_dec=[bits_rx_dec;received(blocklength+1:end)]; %to take only the information bits (not the redudancy)
% %             %here we take only the second part of each decoded block since it
% %             %contains the message without redundancy
% %         end
% 
% 
% 
%         % Check error
%         BER(i,p) = bit_error_rate(bits_rx, bits_with_pilots);
%         %BER(i,2) = bit_error_rate(encoded_message, bits_rx);
%         %BER_moyen(i)=BER_moyen(i)+bit_error_rate(bits, bits_rx);
%         %end
%         %BER_moyen(i)=BER_moyen(i)/5;
        i
    end
end

diff_CFO=CFO-delta_f_tild;
figure;plot(1:length(CFO),diff_CFO);

% figure;
% for p=1:length(CFO)
%     semilogy(10*log10(EbN0),BER(:,p));
%     %figure(25);semilogy(10*log10(EbN0),BER_moyen);
%     xlabel("Eb/N0 (dB)");
%     ylabel("Bit error rate");
%     hold on;
% end