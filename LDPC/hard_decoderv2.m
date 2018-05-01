function out=hard_decoder(H,y,maxit,subH_cell)
%subH_cell is a cell containing the subH matrices, to gain time 

    it=0;
    [row,col]=size(H);
    
    %Create a cell of variables nodes, containing the original value and the
    %number of time a cheker node decides that it should be a one or a zero.
    %c{k} = [valeur,nbre0,nbre1]
    c={};
    for k=1:length(y)
        c{k}=[y(k) 0]; %[valeur,counter]
    end
    
    nextInput = y;
    
    %a changer : 
    %Ne pas cumuler les réponses d'une itération à l'autres !!!
    %plutot faire un compteur -1 +1 -> faire la somme pour la décision
    % Garder en mémoire les réponses des check nodes vers les variables
    % nodes
    % car il faut calculer la réponse des variable nodes vers les check
    % nodes sans tenir compte de la dernière réponse
    
    %nextInputMat = zeros();
    %out_save = [];
    
    while it<maxit
        it = it + 1;
        %out_save = out;
        
        for k=1:length(y)
            c{k}(2)=0; %[valeur,counter]
        end
        
        %for k=1:length(y)
        %    c{k}(1)=nextInput(k); %[valeur,nbre0,nbre1]
        %end
        
        % We do a loop for every Hij which is non null (corresponds to one
        % message that need to be sent back to the validation node j)
        for i=1:row     
            for j=1:col
                if H(i,j)==1
                    
                    %calculate syndrome WITHOUT the validation node j (tip:
                    %we do it only for the check node i to save time.)
                    %Don't forget to discard the validation node j !
                    %subH = horzcat(H(i,1:j-1),H(i,j+1:end));
                    subH=subH_cell{i,j};
                    subY = horzcat(nextInput(1:j-1),nextInput(j+1:end));
                    
                    parcheck = mod(subH*subY',2);
                    
                    if parcheck==0
                       c{j}(2) = c{j}(2) - 1;
                    else                       
                        c{j}(2) = c{j}(2) + 1;
                    end
                end
            end
        end
        %c{:}
        out = zeros(1,length(y));

        %Majority rule
        for k=1:col
            if c{k}(2)<0
                decision=0;
            elseif c{k}(2)>0
                decision=1;
            else
                decision=c{k}(1); %if equal number of 0 or 1, take the received value 
            end
            out(k) = decision;
        end
        
        
%         nextInput = repmat(out,col,1);
%         for e = 1:row
%            nextInput(,e) =  out_save(e);
%         end
        
        nextInput = out;

%         if nnz(mod(H*out',2))==0
%             break
%         end
        if nnz(mod(out*H',2))==0
            break
        end
    end

end