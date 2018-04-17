function out=hard_decoder(H)

[r,c]=size(H);

for check=1:r
    temp=H(check);
    count=0;
    for i=1:length(temp)
        if temp(i)==1
            count=count+1;
        end
        count=mod(count,2);
        if count==1
        
        end
    end
end
            
    
    