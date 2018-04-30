clear all;
close all;

addpath('..')

H = [1 0 0 1 0; 1 1 0 1 0; 0 1 0 0 1; 1 0 1 0 0; 1 1 0 0 0; 0 1 0 1 1; 0 1 0 1 0; 1 0 1 0 1; 0 0 1 1 0; 0 0 1 0 1];
H = H';

Hsub = H(:,1:5);

%useful website
%http://www.di-mgt.com.au/cgi-bin/matrix_stdform.cgi

% H_prim (standard form)
% 1	 0	0	1	1	1	0	0	0	0	
% 1	 1	1	1	0	0	1	0	0	0	
% 1	 1	0	1	1	0	0	1	0	0	
% 1	 0	1	0	1	0	0	0	1	0	
% 0	 1	1	0	0	0	0	0	0	1

H = [1 0 0 0 0 1 0 0 1 1; 0 1 0 0 0 1 1 1 1 0; 0 0 1 0 0 1 1 0 1 1; 0 0 0 1 0 1 0 1 0 1; 0 0 0 0 1 0 1 1 0 0]
%H = [1 0 0 1 1 1 0 0 0 0; 1 1 1 1 0 0 1 0 0 0; 1 1 0 1 1 0 0 1 0 0; 1 0 1 0 1 0 0 0 1 0; 0 1 1 0 0 0 0 0 0 1]
I = H(:,1:5)
Pt = H(:,6:10)

G = [Pt' I]

mod(G*H', 2) %check property

errors=0;


fsymbol=1e6;
T=1/fsymbol;

U=4;
fs = fsymbol*U;
ts = 1/fs;

M = 2;
bits_per_symbol = log2(M)

bits = randi(2,[1,5]) - ones(1,5);
encoded_message = mod(bits*G,2);
original=encoded_message;


%% MAPPING
modulation = 'pam';
%modulation = 'qam';

symb_tx = mapping(encoded_message',bits_per_symbol,modulation);


%% OVERSAMPLING = replicating each symbol U times
U=4;
%figure;
%stem(symb_tx)
symb_tx = upsample(symb_tx,U);



%% Loop for different bit energies +  calculating BER
EbN0 = 10^(0.0001)
    %BER_moyen(i)=0;
    %for j=1:10
    % First half root filter
    beta = 0.3; %imposed
    symb_tx_filtered = halfroot_opti_v2(symb_tx,beta,T,fs);
    %figure;
    %stem(symb_tx_filtered)
    
    
    % Adding noise
    symb_tx_noisy = AWNG(symb_tx_filtered,EbN0,M,fs,modulation,length(encoded_message));
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
    
    real_bits_rx = demapping(symb_tx_noisy,bits_per_symbol,modulation);
    original_errors = nnz(original'-real_bits_rx)

    decoded_bits_rx = soft_decoder_log(H,symb_tx_noisy,20,0.5);

    errors = errors + nnz(decoded_bits_rx-original')