function out=oversample(in)
% in = COLUMN vector

out=zeros(1,2*length(in));

j=0;
for i = 1:length(in)
    j=j+1;
    out(j)=in(i);
    j=j+1;
    out(j)=in(i);
end
out=out';

