function BER = bit_error_rate(bits0,bitsf)
    diff = abs(bits0 - bitsf);
    errors = nnz(diff);
    BER = errors/length(bits0);
end