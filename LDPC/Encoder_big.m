clear all;
close all;

errors = 0;
H0 = makeLdpc(128,256,0,1,3);   %128*256 H -> encodes 128 bits

%% Computing different parts of H with removed columns


for i=i:1000
    bits = randi(2,[128,1]) - ones(128,1);

    [codedbits, H] = makeParityChk(bits, H0, 0);
    
    [row,col]=size(H);
    subH_cell=cell(row,col);
    for i=1:row
        for j=1:col
            if H(i,j)==1
                %calculate syndrome WITHOUT the validation node j (tip:
                %we do it only for the check node i to save time.)
                %Don't forget to discard the validation node j !
                subH_cell{i,j} = horzcat(H(i,1:j-1),H(i,j+1:end));
            end
        end
    end
    
    encoded_message = [codedbits; bits];
    original = encoded_message;

    x = randi(256);
    x2 = randi(256);
    encoded_message(x) = 1-encoded_message(x);
    encoded_message(x2) = 1-encoded_message(x2);

    received = hard_decoderv2(H,encoded_message',50,subH_cell);

    errors = errors + nnz(received'-original);
end

errors