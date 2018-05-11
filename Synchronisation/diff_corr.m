function delta_f_tild=diff_corr(symb_rx,pilot,K,CFO,T)
% attention T is symbol duration
% k is the width of the time difference for the differential correlation
% attention the n in the formula slide 56 is 1,2,3,4,5 and not
% 100,200,300... because this function comes after Gardner so after the
% downsampling

N=length(pilot); %N is the pilot length 
k=1:K;
Dk=zeros(K,N); %each row of the matrix will be for a certain k
for i=1:K
    for n=1:N
        Dk(i,n)=(1/(N-i))*exp(-1j*CFO*i*T).*sum(conj(symb_rx(n+i:n+N-1)).*conj(symb_rx(n:n+N-1-i)).*pilot(1+i:N).*conj(pilot(1:N-i)));
    end
end
%now we have to take the VERTICAL SUM, i.e. the sum on k for a given n
temp1=abs(Dk);
temp1=sum(temp1,2); %vertical sum --> gives a vector
[~,n_tild]=max(temp1);

% temp2=angle(Dk(:,n_tild));
% fact = linspace(2*pi*T,2*pi*T*K,K);
te = 0;
for j=1:K
    te =te + angle(Dk(j,n_tild))/(2*pi*T*j);
end
delta_f_tild = (-1/K)*te;
% temp2 = temp2./(fact.');
% delta_f_tild=(-1/K)*sum(temp2); %=CFO estimate

