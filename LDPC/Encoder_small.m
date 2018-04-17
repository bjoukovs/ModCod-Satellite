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
code = mod(bits*G,2)



