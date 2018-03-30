function symb_tx_noisy = AWNG(symb_tx,EbN0, M, Fsamp, modulation,Nbits)


%% Noise of the channel, AdditiveWhiteNoiseGaussian
% It is a random variable added on each bit during the communication


% INPUTS:
% - symb_tx = symb_txI + j*symb_txQ : vector of ouput symbols (variance 1)
%amplitude Am= (2m - 1 - M) with 
%
% OUTPUTS:
% - symb_tx_noisy : vector of ouput symbols (variance 1)

%% Evaluation of the noise value

Nbpsymb = log2(M);
%Nbits = length(symb_tx)*Nbpsymb; % ATTENTION Nbits doit etre le nbre de
%bits de base, donc en prenant symb_tx qui a été oversamplé ça pose
%problème, donc on le passe en argument avec length(bits) où bits est la
%séquence de bits de base

%Energy = integrale du signal au carré
SignalEnergy = 1/Fsamp*trapz(abs(symb_tx).^2); %this is power at baseband. power at RF must be devided by 2
%Bit energy, *0.5
Eb = 0.5*(SignalEnergy/Nbits);

%Noise power
N0 = Eb / EbN0; %(per hertz)

NoisePower = 2*N0*Fsamp;
% if strcmp(modulation,'pam')
%     NoisePower = NoisePower/2; %Because we only simulate the real part
% end

symb_tx_noisy=zeros(length(symb_tx),1);
%% Addition of the noise
for  i = 1 : length(symb_tx)
    if strcmp(modulation, 'qam')==1
        noiseValue = sqrt(NoisePower/2) * (randn(1, 1) + j*randn(1,1));
    elseif strcmp(modulation, 'pam')==1
        noiseValue = sqrt(NoisePower/2) * (randn(1, 1)); %  METTRE PARTIE IMAGINAIRE QUAND MEME POUR LA SUITE
    end
    symb_tx_noisy(i,1) = symb_tx(i,1) + noiseValue; 
end