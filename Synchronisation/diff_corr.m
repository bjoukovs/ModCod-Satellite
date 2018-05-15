function [delta_f_tild, n_tild]=diff_corr(symb_rx,pilot,K,CFO,T)
% attention T is symbol duration
% k is the width of the time difference for the differential correlation
% attention the n in the formula slide 56 is 1,2,3,4,5 and not
% 100,200,300... because this function comes after Gardner so after the
% downsampling

N=length(pilot); %N is the pilot length 

Dk=zeros(K,N); %each row of the matrix will be for a certain k
for k=0:K-1
    for n=0:N-1
        Dk(k+1,n+1)=(1/(N-k))*exp(-1j*2*pi*CFO*k*T).*sum(conj(symb_rx(n+k+1:n+N)).*symb_rx(n+1:n+N-k).*pilot(k+1:N).*conj(pilot(1:N-k)));
        %Dk(k+1,n+1)=(1/(N-k))*sum(conj(symb_rx(n+k+1:n+N)).*symb_rx(n+1:n+N-k).*pilot(k+1:N).*conj(pilot(1:N-k)));
    end
end
%now we have to take the VERTICAL SUM, i.e. the sum on k for a given n
temp1=abs(Dk);
temp1=sum(temp1,1); %vertical sum --> gives a vector
[~,n_tild]=max(temp1);

% temp2=angle(Dk(:,n_tild));
% fact = linspace(2*pi*T,2*pi*T*K,K);
te = 0;
for j=1:K-1
    te =te + phase(Dk(j+1,n_tild))/j;
end
delta_f_tild = (-1/(K))*te/(2*pi*T);
% temp2 = temp2./(fact.');
% delta_f_tild=(-1/K)*sum(temp2); %=CFO estimate

