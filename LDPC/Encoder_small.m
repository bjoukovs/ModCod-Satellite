clear all;
close all;

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

for i=1:1
    
    
    bits = randi(2,[1,5]) - ones(1,5);
    code = mod(bits*G,2);
    original=code;
    x = randi(10);
    code(x) = 1-code(x);
    %received = hard_decoderv2(H,code,50);
    received = soft_decoder_log(H,code',10,1);
    received = received;
    
    errors = errors + (nnz(mod(received-original,2)) ~= 0);
    
end

errors

%ATTENTION ne marche que si l'erreur a converti un bit 1 en bit 0, pas dans
%l'autre sens --> A CORRIGER



