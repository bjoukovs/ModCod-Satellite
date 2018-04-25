function out=soft_decoder(H,y,maxit,sigma)

    %Translate 0 into -1
    for i=1:length(y)
       
        if y(i)==0
            y(i)=-1;
        end
        
    end

     [row,col]=size(H);
     
     c = y;
     
     it=0;
     
     R0 = zeros(col,row);
     R1 = zeros(col,row);
     Q0 = zeros(row,col);
     Q1 = zeros(row,col);

     PI = zeros(col,1);

     Q_descision0 = zeros(col,1);
     Q_descision1 = zeros(col,1);

     %Step 1
     PI = 1./(1+exp(2*y/sigma^2));

     for i=1:row
        for j=1:col
            if H(i,j)==1

               Q0(i,j) = 1-PI(j);
               Q1(i,j) = PI(j);

            end         
        end
     end
     
     
     % MAIN LOOP
     while it<maxit
         it = it+1

        %Step 2
        for i=1:row
            for j=1:col
                if H(i,j)==1

                    subQ = (ones(1,col-1) - [ Q1(i,1:j-1) Q1(i,j+1:end) ]);
                    subQ = nonzeros(subQ);
                    if isempty(subQ)
                       subQ = [1]; 
                    end

                    R0(j,i) = 0.5 + 0.5*prod(subQ);
                    R1(j,i) = 1 - R0(j,i); 
                end
            end
        end

        %Step 3
        for i=1:row
           for j=1:col
               if H(i,j)==1

                   subR0 = nonzeros(horzcat( R0(j, 1:i-1), R0(j,i+1:end)));
                   subR1 = nonzeros(horzcat( R1(j, 1:i-1), R1(j,i+1:end)));

                   if length(subR0) == 0
                      subR0 = [1];
                   end

                   if length(subR1)==0
                       subR1 = [1];
                   end

                   pi = PI(j);
                   Q0(i,j) = (1-pi)*prod(subR0);
                   Q1(i,j) = pi*prod(subR1);

                   factor = 1/(Q0(i,j)+Q1(i,j));
                   Q0(i,j) = Q0(i,j)*factor;
                   Q1(i,j) = Q1(i,j)*factor;

               end

           end
        end

        %Step 4
        for i=1:col

            subR0 = nonzeros(R0(i,:));
            subR1 = nonzeros(R1(i,:));

            if length(subR0) == 0
              subR0 = [1];
            end

            if length(subR1)==0
               subR1 = [1];
            end
            
            Q_descision0(i) = (1-PI(i))*prod(subR0);
            Q_descision1(i) = PI(i)*prod(subR1);

            factor = Q_descision0(i) + Q_descision1(i);
            Q_descision0(i) = Q_descision0(i)/factor;
            Q_descision1(i) = Q_descision1(i)/factor;

            if Q_descision1(i)>0.5
                c(i) = 1;
                out(i) = 1;
            else
                c(i) = -1;
                out(i) = 0;
            end            

        end
        
        if nnz(mod(H*out',2))==0
            break;
        end
        
     end
     
     %out = c;


end
