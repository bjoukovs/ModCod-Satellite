function dk=diff_corr(symb_rx,pilot,k,CFO,T)
% attention T is symbol duration
% k is the width of the time difference for the differential correlation
% attention the n in the formula slide 56 is 1,2,3,4,5 and not
% 100,200,300... because this function comes after Gardner so after the
% downsampling

N=length(pilot); %N is the pilot length 

dk=zeros(1,N);

for n=1:length(dk)
    Dk(n)=(1/(N-k)).*exp(-1j*CFO*k*T).*sum(conj(symb_rx(n+k:n+N-1)).*conj(symb_rx(n:n+N-1-k)).*pilot(1+k:N).*conj(pilot(1:N-k)))
end