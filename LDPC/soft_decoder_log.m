function out=soft_decoder_log(H,y,maxit,sigma,collocs, collocs_excl, rowlocs, rowlocs_excl)

    %Definition of phi function
    phi = @(x) log10((exp(x)+1)./(exp(x)-1));
    
    %ATTENTION PRENDRE DIRECTEMENT LES SYMBOLES BPSK (partie réelle)
    %(-2*symboles) erreur sur le signe, prendre le moins

     [row,col]=size(H);
     
     c = y;
     
     it=0;
     
     Q = zeros(row,col);
     R = zeros(col,row);
     
     %Initialise step
     for j=1:row %for each check node
         for i=1:col %for each verification node
            if H(j,i)==1
               Q(j,i) = -2*y(i)/sigma^2;
               LCI(i) = Q(j,i);
            end
         end
     end
    
     
     
     %Start iterating
     while it<maxit
        it = it+1;
        
        BETA = abs(Q);
        ALPHA = sign(Q);
        
        %Decision
        LQ = zeros(col,1);
        
        %step 2 : update the check nodes responses
        for j=1:row     %for each check node j
            parfor i=1:length(collocs{j})  %for each verification node connected to the check node j
               
               %sub_collocs = [collocs{j}(1:i-1) collocs{j}(i+1:end)];
               sub_collocs = collocs_excl{j,i};
               sum1 = sum(phi(BETA(j,sub_collocs)));
               
               p1 = prod(nonzeros(ALPHA(j,sub_collocs)));
               
               R(collocs{j}(i),j) = p1*phi(sum1);
                
            end
        end
        
        %step3 : update the verification nodes responses to the check nodes
        for i=1:col     %for each verification node
           for j=1:length(rowlocs{i})   %for each check node connected to the verification node i
              
               %sub_rowlocs = [rowlocs{i}(1:j-1) rowlocs{i}(j+1:end)];
               sub_rowlocs = rowlocs_excl{i,j};
               sum1 = sum(R(i, sub_rowlocs));
               Q(rowlocs{i}(j),i) = LCI(i) + sum1;
               
               
           end  
        end
        
        %step4 : Decision
        for i=1:col
           LQ(i) = LCI(i) + sum(R(i,:));
        end
        
        %step5 : Decision
        for i=1:col
           if LQ(i)<0
              c(i) = 1;
              %out(i)=1;
           else
               c(i)=-1;
               %out(i)=0;
           end
        end
    
         % DEMAPPING
        out = demapping(c',1,'pam');
        
        if nnz(mod(H*out,2))==0
            break;
        end
        
     end
    


end
