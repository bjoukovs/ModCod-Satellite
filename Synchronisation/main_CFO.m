%% Same main as test_mapping_demapping_filter_v2 but with CFO included
%% ATTENTION here no channel coding
clear all;
close all;
addpath('..');

format long

fsymbol=1e6;
T=1/fsymbol;

U=50;
fs = fsymbol*U;
ts = 1/fs;

M = 16;
bits_per_symbol = log2(M);

blocklength=128;
bits = randi(2,bits_per_symbol*50*blocklength,1); %100k symbols
% ATTENTION put more than bits_per_symbol*100*blocklength
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;
bits = bits -1;

%% MAPPING
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% OVERSAMPLING = replicating each symbol U times
symb_tx = upsample(symb_tx,U);

%% Loop for different bit energies +  calculating BER
%EbN0 = logspace(-0.4,2,40);
EbN0 = logspace(0.6,2,10);

%%%%% To investigate CFO only %%%%%
% CFO=0:10:90;
% CFO=deg2rad(CFO);
f_carrier=2e9; %2GHz
CFO = [0 2 10 50]; %in ppm

CFO = f_carrier/1e6*CFO;
phi0=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% To investigate phi0 only %%%%%
% CFO=0;
% phi0=0:10:90;
% phi0=deg2rad(phi0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BER = zeros(length(EbN0),length(CFO));

t=[0:ts:(length(symb_tx)-1)*ts]';

%for p=1:length(phi0)
for p=1:length(CFO)
    p
    for i=1:length(EbN0)
        
        %BER_moyen(i)=0;
        %for j=1:10
        % First half root filter
        beta = 0.3; %imposed
        symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs,U);
        
        % Adding noise
        [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));

        %% Reception: multiplying with exp(j*CFO*t+phi0) to take the effect CFO and phase offset (cf. slide 10)
        t=[0:ts:(length(symb_tx_noisy)-1)*ts]';
        %%%%% Introduce CFO only %%%%%
        symb_tx_noisy = symb_tx_noisy.*exp(2*pi*j.*(CFO(1,p).*t+phi0));
        % 1 symbol is taken at each sampling time ts
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Introduce phi0 only %%%%%
        %symb_tx_noisy = symb_tx_noisy.*exp(1j.*phi0(1,p));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Second half root filter
        beta = 0.3; %imposed

        symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs,U);

        %%%% Remove the extra samples due to the 2 convolutions %%%%
        RRCTaps=25*U+1;
        symb_tx_noisy=symb_tx_noisy(RRCTaps:end-RRCTaps+1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        t = linspace(ts*(RRCTaps/2), (length(symb_tx_noisy)+RRCTaps/2-1)*ts, length(symb_tx_noisy));
        t=t';
        
        %%%% Compensate rotation but not ISI due to modified filter (cf. slide 14) %%%%
        %symb_tx_noisy = symb_tx_noisy.*exp(-1j.*(CFO(1,p).*t(RRCTaps/2:end-3/2*RRCTaps+1)));
        symb_tx_noisy = symb_tx_noisy.*exp(-2*pi*j.*(CFO(1,p).*t));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % DOWNSAMPLING
        symb_tx_noisy = downsample(symb_tx_noisy,U);
        %figure;
        %stem(symb_tx_noisy)

        % DEMAPPING
        if strcmp(modulation,'pam')
            symb_tx_noisy=real(symb_tx_noisy) %if modulation is pam the demapping cannot take complex number
        end
        bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);

        % Check error
        BER(i,p) = bit_error_rate(bits_rx, bits);
        i

    end
end

figure;scatter(real(symb_tx_noisy),imag(symb_tx_noisy));
%figure;plot(symb_tx_noisy)

figure;
for p=1:length(CFO)
    semilogy(10*log10(EbN0),BER(:,p));
    %figure(25);semilogy(10*log10(EbN0),BER_moyen);
    xlabel("Eb/N0 (dB)");
    ylabel("Bit error rate");
    hold on;
end