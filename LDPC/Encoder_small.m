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
I = H(:,1:5)
Pt = H(:,6:10)

G = [Pt' I]

mod(G*H', 2) %check property

bits = randi(2,[1,5]) - ones(1,5)
code = [1     1     1     0     0     0     0     0     1     0]
%mod(bits*G,2)
code = [1     1     1     0     1     0     0     0     1     0]
%mod(code(5)+1,2)
received = hard_decoder(H,code)

%ATTENTION ne marche que si l'erreur a converti un bit 1 en bit 0, pas dans
%l'autre sens --> A CORRIGER



