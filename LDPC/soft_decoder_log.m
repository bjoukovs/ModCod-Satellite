function out=soft_decoder_log(H,y,maxit,sigma)

    %Definition of phi function
    phi = @(x) log10((exp(x)+1)./(exp(x)-1));
    
    %Translate 0 into -1
    for i=1:length(y)

        if y(i)==0
            y(i)=-1;
        end

    end

     [row,col]=size(H);
     
     c = y;
     
     it=0;
     
     Q = zeros(row,col);
     R = zeros(col,row);
     
     %[R] Creating the set of column location of the 1 in H, for each row
     collocs = {}
     for i=1:row
         collocs{i} = [];
        for j=1:col
            if H(i,j)==1
               collocs{i} = [collocs{i} j]; 
            end
        end
     end
     
     %[C] Creating the set of row location of the 1 in H, for each column
     rowlocs = {}
     for j=1:col
         rowlocs{j} = [];
        for i=1:row
            if H(i,j)==1
               rowlocs{j} = [rowlocs{j} i]; 
            end
        end
     end
     
     %Initialise step
     for i=1:row
         for j=1:col
            if H(i,j)==1
               Q(i,j) = 2*y(i)/sigma^2; 
            end
         end
     end
     
     
     %Start iterating
     while it<maxit
        it = it+1;
        
        BETA = abs(Q);
        ALPHA = sign(Q);
        
        %Log-likelyhood ratio
        PI = 1./(1+exp(2*c./sigma^2));
        LCI = log10((1-PI)./PI);
        
        %Decision
        LQ = zeros(col,1);
        
        %step 2
        for i=1:row
            for j=1:length(collocs{i})
               
               sub_collocs = [collocs{i}(1:j-1) collocs{i}(j+1:end)];
               sum1 = sum(phi(BETA(i,sub_collocs)));
               p1 = prod(nonzeros(ALPHA(i,sub_collocs)))
               
               R(collocs{i}(j),i) = p1*phi(sum1);
                
            end
        end
        
        %step3
        for j=1:col
           for i=1:length(rowlocs{j})
              
               sub_rowlocs = [rowlocs{j}(1:i-1) rowlocs{j}(i+1:end)];
               sum1 = sum(R(j,sub_rowlocs)); % TO VERIFY
               Q(rowlocs{j}(i),j) = LCI(i) + sum1;
               
               
           end  
        end
        
        %step4
        for i=1:col
           LQ(i) = LCI(i) + sum(R(i,rowlocs{i})) 
        end
        
        %step5
        for i=1:col
           if LQ(i)<0
              c(i) = 1;
              out(i)=1;
           else
               c(i)=-1;
               out(i)=0;
           end
        end
    
        
        if nnz(mod(H*out',2))==0
            break;
        end
        
     end


end
