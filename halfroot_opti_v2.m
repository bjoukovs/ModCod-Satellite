function out = halfroot_opti_v2(symbup,beta,T,fs)
% fs is the sampling frequncy of the already upsampled signal (symbup)
% RRCTaps = N ? No, it is a paramater that we choose to adjust the
% precision of the filter

RRCTaps = 101; %parameter

if mod(RRCTaps,2)==0
    RRCTaps=RRCTaps+1; %ensure that number of taps is odd
end

%% Frequency Grid
stepOffset=fs/RRCTaps;
highestFreq=stepOffset*(RRCTaps-1)/2;
freqGrid=linspace(-highestFreq,highestFreq,RRCTaps);

%% Construct filter H
H=zeros(length(freqGrid),1); 
% Creation of H by making the correspondance between the frequencies and
% the index of the loop (index i) with freqGrid(i) 
for i=1:RRCTaps
    if abs(freqGrid(i))<=(1-beta)/(2*T)
        H(i,1)=T;
    elseif abs(freqGrid(i))<=(1+beta)/(2*T)
        H(i,1) = (T/2)*(1+cos((pi*T/beta)*(abs(freqGrid(i))-(1-beta)/(2*T))));
    end
end
%the rest is made of zeros

%% Normalizing and taking the square root
h = ifft(H,'symmetric');
h = h/h(1);
H = fft(h);
H = sqrt(H); %halfroot

figure(11); plot(freqGrid,H);

%% Going in temporal domain
Delta_t = 1/fs;
t = linspace((-(RRCTaps-1)/2)*Delta_t,((RRCTaps-1)/2)*Delta_t,RRCTaps);
% shift h of 1 unit to the left with t=[t(1) t] and h=[h;h(end)], because
% small error of shift
%t=[t(1) t];
h = ifft(ifftshift(H),'symmetric');
%h=[h;h(end)];

%% Check if filter is ok by seeing if zeros crossings are at each T
% h_tot=conv(h,h,'same');
% figure(12);plot(t,h_tot); hold on;
% T_check=-T*3:T:T*3;
% plot(T_check,zeros(length(T_check)),'-'); hold on;
% plot(T_check,zeros(length(T_check)),'x');

out = conv(symbup, h, 'same');
    