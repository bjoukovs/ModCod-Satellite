function symb_tx_noisy = AWNG(symb_tx)


%% Noise of the channel, AdditiveWhiteNoiseGaussian
% It is a random variable added on each bit during the communication


% INPUTS:
% - symb_tx = symb_txI + j*symb_txQ : vector of ouput symbols (variance 1)
%amplitude Am= (2m - 1 - M) with 
%
% OUTPUTS:
% - symb_tx_noisy : vector of ouput symbols (variance 1)

%% Evaluation of the noise value
M = 16;
%pam
Nbps = log2(M);
sigma = sqrt(sum(([0:2^Nbps-1]-(2^Nbps-1)/2).^2)/2^Nbps); 
d = 1/sigma;


amp = [];
%d = ;%2d is the distance between 2adj symbol
for m = 1:M-1
   amp(m) =  (2*m - 1 - M)*d;
end
maxAmp = max(amp);
minAmp = min(amp);

diff = maxAmp - minAmp;
noiseAmount = 20;%
maxNoiseAmp = noiseAmount*diff/100;
a = maxNoiseAmp;
b = -maxNoiseAmp;

%% Addition of the noise
for  i = 1 : length(bitStream)
    noiseValue = a + (b-a).*rand(100,1);
    symb_tx_noisy(i) = symb_tx(i) + noiseValue; 
end