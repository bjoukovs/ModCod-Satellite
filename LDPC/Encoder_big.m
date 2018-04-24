clear all;
close all;

errors = 0;
H0 = makeLdpc(128,256,0,1,3);   %128*256 H -> encodes 128 bits

for i=i:1000
    bits = randi(2,[128,1]) - ones(128,1);

    [codedbits, H] = makeParityChk(bits, H0, 0);

    encoded_message = [codedbits; bits];
    original = encoded_message;

    x = randi(256);
    encoded_message(x) = 1-encoded_message(x);

    received = hard_decoderv2(H,encoded_message',50);

    errors = errors + (nnz(received'-original) ~= 0);
end

errors