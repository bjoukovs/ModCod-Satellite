function out_t=halfrootfilter(symbup,cutf,rolloff, fs, M)
%% TO DO:
% use a non rectangular window and take roll-off into account

%% I/Os
% INPUTS:
% symbup = column vector of symbols, already oversampled
% cutf = cutoff frequency [Hz]
% rolloff = roll-off factor
% fs = sampling frequency

% OUTPUTs:
% out_t = column vector of filtered symbols in time domain

%% Rectangular window  (sinc in temporal domain)
% mapping to discrete frequency: f_n = n * fs / N 
figure;stem(symbup);
N = length(symbup);
filter_f = zeros(N,1);
filter_f(1:N/M,1) = 1; %rectangle of amplitude sqrt(ts) that stops at cutoff
%stem(filter_f);

symb_f = fft(symbup,N);
fftshift(symb_f);

%% Multiplication in frequency domain, then inverse Fourier transform
out_f = symb_f .* filter_f;
%out_f = ifftshift(out_f);
out_t = ifft(out_f,N);
figure;stem(out_t);






    