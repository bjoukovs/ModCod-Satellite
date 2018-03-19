function symb_out_t = halfroot_opti(symbup,M,beta,ts)
%% I/Os
% INPUTS:
% symbup = column vector of symbols, already oversampled
% M = oversampling factor
% beta = roll-off
% ts = sampling time

% OUTPUTS:
% symb_out_t = filtered symbols, temporal domain

%% Raised cosine creation in frequency
N = length(symbup);
%transformation factor
fs = 1/ts;
conv = (N/fs); %and not N/fs
filter_f = zeros(N,1);
filter_f(1:floor((1-beta)/(2*ts)*conv),1) = ts;
for i=floor((1-beta)/(2*ts)*conv)+1:floor((1+beta)/(2*ts)*conv)+1  
    filter_f(i,1) = (ts/2)*(1+cos((pi*ts/beta)*(i*fs/N - (1-beta)/(2*ts)))); %ATTENTION i*fs/N and not i*N/fs
end
%filter_f(:,1) = sqrt(filter_f(:,1));

figure(1)
stem(filter_f);

%% Check temporal domain and normalizing
filter_t = ifft(filter_f,'symmetric');
filter_t = filter_t/filter_t(1); %normalizing
filter_f = fft(filter_t);
filter_f = sqrt(filter_f); %halfroot
filter_t = ifftshift(ifft(filter_f,'symmetric'));
figure(2)
stem(filter_t);

%% Filtering
symbup_f = fftshift(fft(symbup));
symb_out_f = symbup_f .* filter_f;
figure(3)
stem(symb_out_f);

symb_out_t = ifftshift(ifft(symb_out_f,'symmetric'));
figure(4)
stem(symb_out_t);
