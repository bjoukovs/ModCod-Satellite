function out=halfroot_opti_v2(symbup,beta,T,fs)
% RRCTaps = N ? No, it is a paramater that we choose to adjust the
% precision of the filter
RRCTaps = 101; %parameter

if mod(RRCTaps,2)==0
    RRCTaps=RRCTaps+1; %ensure that number of taps is odd
end

stepOffset=fs/RRCTaps;
highestFreq=stepOffset*(RRCTaps-1)/2;
freqGrid=linspace(-highestFreq,highestFreq,RRCTaps);

H=zeros(length(freqGrid),1); %we just take the positive frequencies
% H(1:(1-beta)/(2*T),1)=T;
% for f=(1-beta)/(2*T):(1+beta)/(2*T)
%     H(f,1) = (T/2)*(1+cos((pi*T/beta)*(f-(1-beta)/(2*T))));
% end
for i=1:RRCTaps
    if abs(freqGrid(i))<=(1-beta)/(2*T)
        H(i,1)=T;
    elseif abs(freqGrid(i))<=(1+beta)/(2*T)
        H(i,1) = (T/2)*(1+cos((pi*T/beta)*(abs(freqGrid(i))-(1-beta)/(2*T))));
    end
end
%the rest is made of zeros

h = ifft(H,'symmetric');
h = h/h(1);
H = fft(h);
H = sqrt(H); %halfroot
h = ifftshift(ifft(H,'symmetric'));

%figure(10); stem(h);
%figure(10); stem(h);
figure(11); stem(H);

h_tot=conv(h,h);
figure(12);plot(h_tot)

    