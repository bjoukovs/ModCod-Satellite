clear; close all;clc;

fsymbol=1e6;
T=1/fsymbol;

U=1;
fs = fsymbol*U;
ts = 1/fs;
M = 16;
bits_per_symbol = log2(M)

SIZE_PILOT=18; %4*bits_per_symbol; %nbre of symbols in pilot
LENGTH_FRAME=18*2;

bits = randi(2,bits_per_symbol*10*LENGTH_FRAME,1); %100k symbols

bits = bits -1;

%% MAPPING
modulation = 'pam';
%modulation = 'qam';
symb_tx = mapping(bits,bits_per_symbol,modulation);

%% INSERTING PILOTS
pilot=makePilot(modulation,bits_per_symbol,SIZE_PILOT);
symb_tx=intoFrames(symb_tx,LENGTH_FRAME,pilot); %attention now symb_tx is longer
bits_with_pilots=demapping(symb_tx,bits_per_symbol,modulation); %save the bits of the symbol containing pilots,

EbN0 = 100;
%%%%% To investigate CFO only %%%%%
f_carrier=2e9; %2GHz
CFO=0:10^(-6)*f_carrier:10*10^(-6)*f_carrier;
phi0=0;
delta_f_hat=zeros(1,length(CFO));
Tp = 1/(2*f_carrier);
for p=1:length(CFO)
    p
    t=[0:ts:(length(symb_tx)-1)*ts]';
    symb_tx = symb_tx.*exp(1j.*(CFO(1,p).*t+phi0));

    check_length=SIZE_PILOT+LENGTH_FRAME;
    K=8;
    b1=1;
    b2=check_length+1;
    delta_f_hat(1,p) = diff_corr(symb_tx(b1:b2),pilot,K,CFO(1,p),Tp);
end  
diff_CFO=CFO-delta_f_hat;
figure;plot(1:length(CFO),diff_CFO);
figure;stem(1:length(CFO),delta_f_hat);
