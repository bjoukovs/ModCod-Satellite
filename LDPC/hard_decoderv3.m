function out=hard_decoderv2(H,y,maxit,subH_cell)
%subH_cell is a cell containing the subH matrices, to gain time 

    it=0;
    [row,col]=size(H);
    
    %Create a cell of variables nodes, containing the original value and the
    %number of time a cheker node decides that it should be a one or a zero.
    %c{k} = [valeur,nbre0,nbre1]
    
    
    c=[];
    counter=[];
    for k=1:length(y)
        if y(k)==0
            c(k) = -1;
            counter(k) = 0;
        else
            c(k) = 1;
            counter(k) = 0;
        end
    end
    
    %R(i,j) = Previous messages from check node j to verification node i
    R = zeros(col, row);
    R_new = R;
    
    
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
        
        counter = zeros(1,col);
        
       for j=1:row %for each check node j
           
           %First, the previous message from check node j to all vnodes is
           %substracted, so the majority rule can be applied
           
           if it == 1
               %Shortcut at first iteration
               vnode_messages = y;
           else
               
               %Normal case -> majority rule without the last R message
               %from check node j
               vnode_messages = c+sum(R,2)'-R(:,j)';

               for k=1:length(vnode_messages)
                  if vnode_messages(k)<0
                     vnode_messages(k) = 0; 
                  elseif vnode_messages(k)>0
                      vnode_messages(k) = 1;
                  else
                      vnode_messages(k) = y(k);
                  end
               end
           end

            for i=1:col %for each verification node i
                if H(j,i)==1 %if i is connected to j
                    
                    %update response of check node j to verification node i
                    %(don't take into account the message from vnode i and all the preivous messages from cnode j to vnode c(i' != i)
                    
                    
                    %To do this, we take all the v nodes connected to j,
                    %except i
                    subH = subH_cell{j,i};
                    subVNodesMessage = [vnode_messages(1:i-1) vnode_messages(i+1:end)];
                    
                    local_syndrome = mod(subH*subVNodesMessage',2);
                    
                    if local_syndrome==0
                        %The vnode should be 0
                        
                        counter(i) = counter(i) - 1;
                        R_new(i,j) = -1;
                        
                        %NEW R !
                        
                    else
                        %The vnode should be 1
                        
                        counter(i) = counter(i) + 1;
                        R_new(i,j) = 1;
                        
                    end                    
                end
            end           
       end
       
       R = R_new;
        
        %Now, for the decision, we apply the majority rule
        for i=1:col
           decision(i) = c(i) + counter(i);
           
           if decision(i) > 0
               %decision(i) = 1;
               out(i) = 1;
           elseif decision(i) < 0
               %decision(i) = -1;
               out(i) = 0;
           else
               %decision(i) = c(i);
               out(i) = y(i);
           end
            
            
        end
        
        if nnz(mod(out*H',2))==0
            break
        end
    end

end

            
        
            
    
 