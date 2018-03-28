function out=convol_lab1ex3(h,x)
N=length(h);
M=length(x);
out=zeros(1,N+M-1); %N+M-1
%h=[zeros(1,N+M-1-N) h zeros(1,N+M-1-N)]; --> pas besoin de padder pour h
%x=[zeros(1,N+M-1-M) x zeros(1,N+M-1-M)]; %padding zeros
%h(N+1:N+M-1,:)=0; 
%x(M+1:N+M-1,:)=0;
    
%L=max(N,M);
% for t=1:length(out) %N+M-1
%     for tau=1:N
%         if t-tau>0
%             out(t)=out(t)+h(tau)*x(t-tau);
%         end
%     end
% end
h=flip(h);
x=[zeros(1,N+M-1-M) x zeros(1,N+M-1-M)];

for t=1:length(x)-(N-1)
    for i=1:length(h)
        out(t)=out(t)+x(i+t-1)*h(i);
    end
end
