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

SIZE_PILOT=4*bits_per_symbol; %nbre of symbols in pilot
LENGTH_FRAME=300-SIZE_PILOT;

bits = randi(2,bits_per_symbol*10*LENGTH_FRAME,1); %100k symbols
bits = bits -1;

%% MAPPING
%modulation = 'pam';
modulation = 'qam';

symb_tx = mapping(bits,bits_per_symbol,modulation);

%% INSERTING PILOTS

pilot=makePilot(modulation,bits_per_symbol,SIZE_PILOT);
symb_tx=intoFrames(symb_tx,LENGTH_FRAME,pilot); %attention now symb_tx is longer

bits_with_pilots=demapping(symb_tx,bits_per_symbol,modulation); %save the bits of the symbol containing pilots,

%% OVERSAMPLING = replicating each symbol U times
symb_tx = upsample(symb_tx,U);

%% Loop for different bit energies +  calculating BER
EbN0 = logspace(0,1.6,8);
CFO = 0;
K = [2 8 16];

for j=1:length(K)
    j
    std_of_error = zeros(1,length(EbN0));
    error = zeros(1,5);
    for i=1:length(EbN0)
        i
        for exp=1:5
            exp
            %1st halfroot filter
            beta = 0.3;
            symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs,U);
            %Noise
            [symb_tx_noisy sigma] = AWNG(symb_tx_filtered,EbN0(1,i),M,fs,modulation,length(bits));
            %2nd half root filter
            symb_tx_noisy = halfroot_opti_v2(symb_tx_noisy,beta,T,fs,U);
            %%%% Removing the extra samples due to the 2 convolutions %%%%
            RRCTaps=25*U+1;
            symb_tx_noisy=symb_tx_noisy(RRCTaps:end-RRCTaps+1);
            % DOWNSAMPLING
            symb_tx_noisy = downsample(symb_tx_noisy,U);

            %Check cor
            check_length=SIZE_PILOT+LENGTH_FRAME;        
            b1=1;
            b2=check_length+1;
            delta_f_tild = diff_corr(symb_tx_noisy(b1:b2),pilot,K(j),CFO,T);
            %No CFO, error = deta
            error(exp,i) = delta_f_tild;
        end
        std_of_error(i) = std(error(:,i));
    end
    figure(1);plot(10*log10(EbN0),std_of_error);title('CFO error as a function of K');
    ylabel('frequency error stdev [ppm]'); xlabel('Eb/N0');hold on;
    legend('K = 2','K = 8','K = 16');
end

