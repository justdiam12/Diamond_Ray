%Hamming
 
clear
N_e=8; Nfft=1024;
d=0.5; f=1500; c=1500;
lambda=c/f;
L=N_e*d; 
beam_width=(lambda/(2*L))*180/pi
m=0:Nfft-1;
x=zeros(1,Nfft);
a=ones(1,N_e)/N_e;
x(1:N_e)=a;
y=fft(x,Nfft);
w1=y(1:Nfft/2);w2=y(Nfft/2+1:Nfft);
w2=fliplr(w2);w1=fliplr(w1);
B=[w1 w2];
B_dB=20*log10(abs(B));
 
figure(1)
theta=real(asin((m-Nfft/2)*lambda/(d*Nfft)));
theta_deg=theta*180/pi;
plot(theta_deg, B_dB);grid;
xlabel('Angle - deg');xlabel('Amplitude - dB)')
axis([ -90 90 -40 5]) 
figure(2);
polar(theta+pi/2,B_dB+40)
