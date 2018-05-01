function [collocs collocs_excl rowlocs rowlocs_excl] = get_row_col_locs(H)

    [row,col]=size(H);
    
    %[R] Creating the set of column location of the 1 in H, for each row
     %collocs_excl{i,j} is the columns locations of 1 in the ith row of H,
     %without the jth element
     collocs = {};
     collocs_excl = {};
     for i=1:row
         temp_collocs = [];
         %collocs{i} = [];
        for j=1:col
            if H(i,j)==1
               temp_collocs = [temp_collocs j]; 
            end
        end
        
        collocs{i} = temp_collocs;
        
        for k=1:length(temp_collocs)
           collocs_excl{i,k} = [temp_collocs(1:k-1) temp_collocs(k+1:end)]; 
        end
        
     end
     
     %Time saving : already get all the sub-collocs
     
     
     %[C] Creating the set of row location of the 1 in H, for each column
     rowlocs = {};
     rowlocs_excl = {};
     for j=1:col
         temp_rowlocs = [];
        for i=1:row
            if H(i,j)==1
               temp_rowlocs = [temp_rowlocs i]; 
            end
        end
        
        rowlocs{j} = temp_rowlocs;
        
        for k=1:length(temp_rowlocs)
           rowlocs_excl{j,k} = [temp_rowlocs(1:k-1) temp_rowlocs(k+1:end)]; 
        end
        
     end

end