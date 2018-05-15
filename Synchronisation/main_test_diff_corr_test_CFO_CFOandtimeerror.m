%% Same main as test_mapping_demapping_filter_v2 but with CFO included
%% ATTENTION here no channel coding

clear all;
close all;
clc;
addpath('..'); %parent directory

format long

fsymbol=1e6;
T=1/fsymbol;

U=100;
fs = fsymbol*U;
ts = 1/fs;

M = 4;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*2800,1); 
%bits = randi(2,bits_per_symbol*10*LENGTH_FRAME,1); %100k symbols
% ATTENTION put more than bits_per_symbol*100*blocklength
%bits = ones(bits_per_symbol*10,1); %1k symbols
%bits(5) = 0;

bits = bits -1;

%% MAPPING
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%size_pilot=[10 20 40];
size_pilot=20;

fcarrier=2e9;
CFO=[0 10*10^(-6)*fcarrier];
sample_time_shift=[0 0.02];

EbN0 = logspace(0,1.6,8);
delta_f_tild_final=zeros(length(CFO),length(EbN0));

NBRE_EXP=5;

delta_f_tild=zeros(NBRE_EXP,length(EbN0));

symb_tx_save=symb_tx;

for n=1:length(CFO)
    
    for p=1:NBRE_EXP
        %SIZE_PILOT=4*bits_per_symbol; %nbre of symbols in pilot
        LENGTH_FRAME=300-size_pilot;
        %% INSERTING PILOTS

        symb_tx=symb_tx_save;
        pilot=makePilot(modulation,bits_per_symbol,size_pilot);
        symb_tx=intoFrames(symb_tx,LENGTH_FRAME,pilot); %attention now symb_tx is longer

        %bits_with_pilots=demapping(symb_tx,bits_per_symbol,modulation); %save the bits of the symbol containing pilots,
        % if comparison is needed after

        %% OVERSAMPLING = replicating each symbol U times
        symb_tx = upsample(symb_tx,U);

        %% Loop for different bit energies +  calculating BER
        for i=1:length(EbN0)
            beta = 0.3; %imposed
            symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs,U);
            % Adding noise
            [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));

            %% Reception: multiplying with exp(j*CFO*t+phi0) to take the effect CFO and phase offset (cf. slide 10)
            %%%%% To investigate CFO only %%%%%
            RRCTaps=25*U+1;
            t=[-(RRCTaps/2)*ts : ts : ((length(symb_tx_noisy)-1)-RRCTaps/2)*ts]';
            symb_tx_noisy = symb_tx_noisy.*exp(1j.*CFO(1,n).*t);
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
             % Second half root filter
             beta = 0.3; %imposed 
             symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs,U);
             %%%% Removing the extra samples due to the 2 convolutions %%%%
             symb_tx_noisy=symb_tx_noisy(RRCTaps:end-RRCTaps+1);
             
            %%%% Adding t0 shift %%%%
            symb_tx_noisy = circshift(symb_tx_noisy,round(sample_time_shift(n)*U));
            %%%%%%%%%%%%%%%%%%%%%%%%%
             
             % DOWNSAMPLING
              symb_tx_noisy = downsample(symb_tx_noisy,U);

            %%%% CFO estimation: we estimate the CFO on one frame %%%%
            check_length=size_pilot+LENGTH_FRAME;
            K=8;
            b1=1;
            b2=check_length+1;
            %delta_f_tild(n,p) = diff_corr(symb_tx_noisy(b1:b2),pilot,K,CFO(1,p),T);
            %delta_f_tild(n,i) = delta_f_tild(n,i) + diff_corr(symb_tx_noisy(b1:b2),pilot,K,0,T);
            delta_f_tild(p,i)=diff_corr(symb_tx_noisy(b1:b2),pilot,K,CFO(1,n),T);

            i
        end
    end
    delta_f_tild_final(n,:)=std(delta_f_tild);
end

delta_f_tild_final_ppm=delta_f_tild_final*1e6/(1e9)

figure;
% for q=1:length(size_pilot)
%     %diff_CFO=abs(zeros(length(delta_f_tild(q,:)))-delta_f_tild(q,:));
%     plot(EbN0,abs(diff_CFO)); hold on;
% end
plot(10*log10(EbN0),abs(delta_f_tild_final_ppm));
xlabel('EbN0 [dB]');
ylabel('Error on CFO');
