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

M = 4;
bits_per_symbol = log2(M)

%SIZE_PILOT = [4 8 16 32];
SIZE_PILOT = 20;

K = 8

%SIZE_PILOT=4*bits_per_symbol; %nbre of symbols in pilot
LENGTH_FRAME=300-SIZE_PILOT;


%% MAPPING
%modulation = 'pam';
modulation = 'qam';


%% EBN0

EbN0 = logspace(-0.1,1,10);
%EbN0 = 100;
%EbN0 = [10^1.5]

f_carrier=2e9; %2GHz
CFO = [0 10];
CFO=CFO*f_carrier/1e6;

epsilon = [0 0.02];

        

n_tild_stdev = zeros(length(EbN0),length(CFO));
cfo_stdev = zeros(length(EbN0),length(CFO));
%n_tild_stdev = zeros(length(EbN0),length(K));

for pl=1:length(SIZE_PILOT)
    for kl=1:length(K)
    


        %%%%% To investigate CFO only %%%%%

        phi0=0;
        delta_f_tild=zeros(1,length(CFO));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        BER = zeros(length(EbN0),length(CFO));

        for p=1:length(CFO)
        %    p

            for i=1:length(EbN0)

                NEXP = 50;
                for xp=1:NEXP

                    bits = randi(2,bits_per_symbol*2*LENGTH_FRAME(pl),1); %100k symbols
                    bits = bits -1;

                    %MAPPING
                    symb_tx = mapping(bits,bits_per_symbol,modulation);

                    pilot=makePilot(modulation,bits_per_symbol,SIZE_PILOT(pl));
                    symb_tx=intoFrames(symb_tx,LENGTH_FRAME(pl),pilot); %attention now symb_tx is longer

                    bits_with_pilots=demapping(symb_tx,bits_per_symbol,modulation); %save the bits of the symbol containing pilots,
                    % if comparison is needed after


                    %% OVERSAMPLING = replicating each symbol U times

                    symb_tx = upsample(symb_tx,U);

                    % First half root filter
                    beta = 0.3; %imposed
                    symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs,U);
                    

                    % Adding noise
                    [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));


                    %% Reception: multiplying with exp(j*CFO*t+phi0) to take the effect CFO and phase offset (cf. slide 10)
                    %%%%% To investigate CFO only %%%%%

                    RRCTaps=25*U+1;
                    t=[-(RRCTaps/2)*ts : ts : ((length(symb_tx_noisy)-1)-RRCTaps/2)*ts]';
                    symb_tx_noisy = symb_tx_noisy.*exp(2*pi*j.*(CFO(1,p).*t+phi0));
                    
                    %Time shift error
                    symb_tx_noisy = circshift(symb_tx_noisy,round(epsilon(p)*U));


                    %%%%%%%%%%%%%%%%%%%%%%%%%%

                     % Second half root filter
                     beta = 0.3; %imposed 
                     symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs,U);
                     %%%% Removing the extra samples due to the 2 convolutions %%%%
                     symb_tx_noisy=symb_tx_noisy(RRCTaps:end-RRCTaps+1);

                     % DOWNSAMPLING
                      symb_tx_noisy = downsample(symb_tx_noisy,U);


                    %%%% CFO estimation: we estimate the CFO on one frame %%%%
                    %%%% There is no need to compensate it because it is too
                    %%%% complicated
                    check_length=SIZE_PILOT+LENGTH_FRAME;
                    %K=5;
                    b1=1;
                    b2=check_length+1;
                    %delta_f_tild(1,p) = diff_corr(symb_tx_noisy(b1:b2),pilot,K,CFO(1,p),T);
                    %[deltaCFO, ntild] = diff_corr(symb_tx_noisy(b1:b2),pilot,K,CFO(1,p),T);
                    [deltaCFO, ntild] = diff_corr(symb_tx_noisy(b1:b2),pilot,K(kl),0,T);

                    n0 = 1; %Expected ToA sample
                    n_tild_stdev(i,p) =  n_tild_stdev(i,p) + (ntild - n0)^2;
                    cfo_stdev(i,p) =  cfo_stdev(i,p) + (deltaCFO - CFO(1,p))^2;

                    xp
                end
                n_tild_stdev(i,p) = sqrt(n_tild_stdev(i,p)/NEXP);
                cfo_stdev(i,p) = sqrt(cfo_stdev(i,p)/NEXP)*f_carrier/1e6;
                i
            end
            
        end


    end
    
end

figure;plot(10*log10(EbN0),n_tild_stdev);
figure;plot(10*log10(EbN0),cfo_stdev);