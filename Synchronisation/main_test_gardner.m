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

U=100;
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
symb_tx_upsampled = upsample(symb_tx,U);
n_original = 1:U:length(symb_tx_upsampled);
sample_time_shift = linspace(0,0.1,5); %expressed in the form of a percentage of the original sampling time

%% Loop for different bit energies +  calculating BER
%EbN0 = logspace(0,2,10);
%EbN0 = [10^10]
EbN0 = [10^2]
%EbN0 = [10^0.5]
BER = zeros(length(EbN0),length(sample_time_shift));
std_vect = zeros(1,length(sample_time_shift));

for p=1:length(sample_time_shift)
    
    n_sampling = n_original + ones(1,length(n_original))*U*sample_time_shift(p);
    for i=1:length(EbN0)
        
        % First half root filter
        beta = 0.3; %imposed
        symb_tx_filtered = halfroot_opti_v2(symb_tx_upsampled,beta,T,fs,U);

        % Adding noise
        [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));
        
%         %%%%% To add the CFO impqct
%         f_carrier=2e9; %2GHz
%         CFO=0:10^(-6)*f_carrier:10*10^(-6)*f_carrier;
%         phi0=0;
%         delta_f_tild=zeros(1,length(CFO));
%         RRCTaps=25*U+1;
%         t=[-(RRCTaps/2)*ts : ts : ((length(symb_tx_noisy)-1)-RRCTaps/2)*ts]';
%         symb_tx_noisy = symb_tx_noisy.*exp(1j.*(CFO(1,p).*t+phi0));
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Second half root filter
        beta = 0.3; %imposed
        symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs,U);
        
        %%%% Removing the extra samples due to the 2 convolutions %%%%
        RRCTaps=25*U+1;
        symb_tx_noisy=symb_tx_noisy(RRCTaps:end-RRCTaps+1);
         
        %%%% Adding t0 shift %%%%
        symb_tx_noisy = circshift(symb_tx_noisy,round(sample_time_shift(p)*U));
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        % DOWNSAMPLING TAKING WITH GARDNER
        kappa = 0.002/2;
        [sampled_signal,eps_ev] = gardner(kappa,symb_tx_noisy,U,1,n_original);
        
        % DEMAPPING
        bits_rx = demapping(sampled_signal.',bits_per_symbol,modulation);

        % Check Gardner: evolution of eps_tild
        figure(5);plot(eps_ev);  title('Convergence of the epsilon for different value of error weight, EbN0 = 5dB');
        ylabel('compensate error weight (sample time shift, % of the symbol frequency)'); xlabel('Samples');hold on;
        legend('eps = 0','eps = 0.025','eps = 0.05','eps = 0.075','eps = 1');
        
        % Standart deviation
        std_vect(p) = std(eps_ev);
        
        % Check error
        BER(i,p) = bit_error_rate(bits_rx, bits);
        i
    end
end

figure; plot(sample_time_shift,std_vect); title('Standard deviation Eb/N0 = 5dB');
xlabel('error weight (sample time shift, % of the symbol frequency)'); ylabel('Standard deviation');
% for p=1:length(sample_time_shift)
%     semilogy(10*log10(EbN0),BER(:,p));
%     %figure(25);semilogy(10*log10(EbN0),BER_moyen);
%     xlabel('Eb/N0 (dB)');
%     ylabel('Bit error rate');
%     hold on;
% end