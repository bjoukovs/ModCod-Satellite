clear all;
close all;

M = 16;
bits_per_symbol = log2(M)

bits = randi(2,bits_per_symbol*100,1);
bits = bits -1;

%% MAPPING
modulation = "pam";
%modulation = "qam";

symb_tx = mapping(bits,bits_per_symbol,modulation)


%% DEMAPPING