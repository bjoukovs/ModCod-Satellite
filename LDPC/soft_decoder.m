function out=soft_decoder(H,y,maxit,sigma)


     [row,col]=size(H);
     
     c = y; 
     
     Q0 = zeros(row,col);
     Q1 = zeros(row,col);
     R0 = zeros(col,row);
     R1 = zeros(col,row);
     
     PI = zeros(row,1);
     
     Q_descision0 = zeros(row,1);
     Q_descision1 = zeros(row,1);
     
     %Step 1
     for i=1:row
        for j=1:col
            if H(i,j)==1

               Q0(i,j) = 1/(1+exp(-2*y(i))/sigma^2);
               Q1(i,j) = 1/(1+exp(2*y(i))/sigma^2);
               
               PI(i) = Q1(i,j);
                   
            end         
        end
        
        %Step 2
        for i=1:row
            for j=1:col
                if H(i,j)==1
                    
                    subQ = (ones(row-1,1) - vertcat( Q1(1:i-1,j), Q1(i+1:end,j) )); %On prend la j colonne, on retire la ligne i
                    
                    R0(j,i) = 0.5 + 0.5*prod(subQ);
                    R1(j,i) = 1 - R0(j,i); 
                end
            end
        end
        
        %Step 3
        for i=1:row
           for j=1:col
               if H(i,j)==1
              
                   subR0 = vertcat( R0(1:j-1, i), R0(j+1:end, i));
                   subR1 = vertcat( R1(1:j-1, i), R1(j+1:end, i));

                   pi = Q1(i,j);
                   Q0(i,j) = (1-pi)*prod(subR0);
                   Q1(i,j) = pi*prod(subR1);

                   factor = 1/(Q0(i,j)+Q1(i,j));
                   Q0(i,j) = Q0(i,j)*factor;
                   Q1(i,j) = Q1(i,j)*factor;
               
               end
               
           end
        end
        
        %Step 4
        for i=1:row
            
            subR = R0(1,:);
            Q_descision0(i) = (1-PI(i))*prod(subR);
            Q_descision1(i) = PI(i)*prod(subR);
            
            factor = Q_descision0(i) + Q_descision1(i);
            Q_descision0(i) = Q_descision0(i)/factor;
            Q_descision1(i) = Q_descision1(i)/factor;
            
            if Q_descision1(i)>0.5
                c(i) = 1;
            else
                c(i)=0;
            end            
            
        end
        
     end
     
     out = c;


end
