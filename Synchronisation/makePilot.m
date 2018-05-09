function output = makePilot(modulation,bits_per_symbol,pilot_length)
    
    %Generates a random pilot of length pilot_length (number of symbols of
    %the pilot) depending on the chosen modulation and the number of bits
    %per symbols
    
    addpath('..');
    
    N = pilot_length*bits_per_symbol;
    bits = randi([0 1],N,1);
    
    output  = mapping(bits,bits_per_symbol,modulation);

    

end