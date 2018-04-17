function out=hard_decoder(H,y)

[row,col]=size(H);

%Create a cell of variables nodes, containing the original value and the
%number of time a cheker node decides that it should be a one or a zero.
%c{k} = [valeur,nbre0,nbre1]
c={};
for k=1:length(y)
    c{k}=[y(k) 0 0]; %[valeur,nbre0,nbre1]
end

%f is a column vector containing the modulo sum of the variable nodes
%linked to the cheker node f(i) (slide 10 modulo sum c0 + c1 + c2 + c5 =? 0)
f=mod(H*y',2)

%This loop implements what is done in slide 13
for i=1:length(f)
    %If f(i) = 1 it means that the modulo sum is different from 0 and then
    %that the bit should change its value acording to this cheker done:
    %increments the nbre0 or nbr1 parameter of c{j} depending on the value
    %of c{j}(1), the original value, the y(j)
    if f(i)==1
        for j=1:col
            if H(i,j)==1
                %to know if we have to increment nbre0 or nbre1, i.e. if we
                %have to increment at index 2 (-->nbre0) or index 3(-->nbre1)
                %of c{j}, we take the index equal to 2 + the modulo of [ the original value c{j}(1)
                %(i.e. y(j)) + 1 ]
                c{j}(2+mod(c{j}(1)+1,2)) = c{j}(2+mod(c{j}(1)+1,2)) + 1;
            end
        end
    else %if f(i)=0 then increment nbre0 if y(j)=0 and nbre1 if y(j)=1, so it means that the checker node
         %tells the variable node to keep its original value
        for j=1:length(H(i))
            if H(i,j)==1
                c{j}(2+mod(c{j}(1),2)) = c{j}(2+mod(c{j}(1),2)) + 1;
            end
        end
    end
end

out = zeros(1,length(y));

%Majority rule
for k=1:length(c)
    decision = round( ( c{k}(1) + c{k}(3) )/ (1 + c{k}(2) + c{k}(3)) ); %moyenne du nombre de 0 et de 1, arrondi
    %if nbre1 > nbre0 then decision = 1 thanks to the rounding, and
    %it means that the node was told more time to be a 1 and so final value
    %of the variable node is a 1. If nbre1 < nbre0 then decision = 0 and
    %final value of variable node is 0.
    out(k) = decision;
end

end

            
        
            
    
 