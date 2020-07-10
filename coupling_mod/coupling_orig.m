%% load data
load('LFP_data.mat')
addpath('/Users/rcruces/matlab_toolbox/PAC')

%% Filter frequencies
Fs=5000;                             %CAlculo de la frecuencia de sampleo
NyFq=Fs/2;                           %CAlulo del ratio de Nyquist
Cut20=20/NyFq;                       %L?mite inferior de gamma
Cut40=40/NyFq;                       %L?mite superior de gamma
Cut4=4/NyFq;                         %L?mite inferior de theta
Cut11=11/NyFq;                       %L?mite superior de theta
[Bg,Ag]=butter(4,[Cut20 Cut40]);     %Coeficientes Butterworth para gamma    
[Bt,At]=butter(2,[Cut4 Cut11]);      %Coeficientes Butterworth para theta
gamma=filtfilt(Bg,Ag,LFP);           %Filtro sin delay gamma
theta=filtfilt(Bt,At,LFP);           %Filtro sin delay theta
 
filteredgamma=dataset(Time,gamma);
filteredtheta=dataset(Time,theta);

% Gamma peaks over theta cycles
plot(Time,gamma,Time,theta)
axis([0 1 -0.3 0.3]) %de 0 a 1 s y de -0.2 a 0.2 mV
shg
clf; clf; 

%% Gamma Hilbert > rotation > peak detection
hilGamma=hilbert(gamma); %transformada de hilbert de la se?al de campo
angGamma=angle(hilGamma); %extrae la fase instant?nea de la se?al de campo en radianes
angGamma=angGamma-2*pi*floor(angGamma/(2*pi));%expande los ?ngulos obtenidos a un CArculo completo (2*pi)
 
negGam=-angGamma;%invierte la se?al para encontrar los 90? = 0 rad
[~,subGamma]=findpeaks(negGam,'MinPeakHeight',-0.5); %encuentra la posici?n todos los picos que superan el threshold 
disp('Total de picos de gamma:')
TotalGPeaks=length(subGamma);

% plot of the peaks
plot(Time, negGam)
hold on
plot(Time(subGamma),negGam(subGamma), 'ro') 
close

%% Theta Hilbert > rotation > peak detection
hilTheta=hilbert(theta); %transformada de hilbert de la se?al de campo
angTheta=angle(hilTheta); %extrae la fase instant?nea de la se?al de campo en radianes
angTheta=angTheta-2*pi*floor(angTheta/(2*pi));%expande los ?ngulos obtenidos a un CArculo completo (2*pi)
 
negTheta=-angTheta;%invierte la se?al para encontrar los 90? = 0 rad
[~,subTheta]=findpeaks(negTheta,'MinPeakHeight',-0.01); %encuentra la posici?n todos los picos que superan el threshold 
disp('Total de picos de Theta:')
TotalTPeaks=length(subTheta);

% PLot of Theta hilbert
plot(hilTheta);

% plot of the peaks
plot(Time, negTheta)
hold on
plot(Time(subTheta),negTheta(subTheta), 'ro') 

%% Total Gamma peaks per theta cycle
disp('Gamma peaks per Theta cycle')
GpeakPerTheta=TotalGPeaks/TotalTPeaks;
 
%encuentra el angulo de theta donde se encuentra los picos de gamma
gammaTheta=angTheta(subGamma);

% CALCULO DE LA LONGITUD DEL VECTOR
disp('Vector unitario medio:')
u=sum(exp(1i*gammaTheta));%calcula el vector unitario medio (f?rmula de Euler)
disp('Longitud total del vector Normalizado:')
absVecLength=abs(u)/TotalGPeaks; %obtiene el radio resultante (r) normalizado entre 0 y 1

%CALCULO DEL ANGULO DEL VECTOR GAMMA OVER THETA (gammaTheta)
disp('?ngulo del vector en pi valores:')
AngGammaTheta=angle(u) % Calcula el angulo del vector normalizado
disp('?ngulo del vector en radianes:')
angGammaTheta_Rad=AngGammaTheta-2*pi*floor(AngGammaTheta/(2*pi)) %Expande el angR a 2pi


figure(1)
subplot(1,2,1)
rose(gammaTheta)
subplot(1,2,2)
compass(u/TotalGPeaks)
 
% figure(2)
% subplot(1,2,1)
% rose(RgamTh)
% subplot(1,2,2)
% compass(sum(exp(1i*RgamTh))/length(RgamTh))
 
clf; clf;

%% Gamma Envelope
% 
%CALCULO DE LA envDemean
[envPeak,~]=envelope(gamma,100,'peak');%CAlculo de la envDemean ajustada al pico con 100 puntos de resoluci?n
envDemean=envPeak-mean(envPeak);%ajuste de la envDemean a la basal
% grafico del envolvente sobre Gamma
plot(Time, gamma, Time, envPeak)

%% TRANSFORMADA DE FOURIER DE AMBAS OSCILACIONES
% Nombra LFP a la variable que contiene el registro de campo
% Nombra Time a la variable que contiene el tiempo de registro
 
%Frecuencia de sampleo (=n? de puntos/tiempo de registro (s))
Fs=5000;
%CAlulo del tama?o de la variable
L=length(envDemean);
%optimizaci?n del tama?o de las variables para la FFT (el m?tiplo de 2 mayor m?s cercano)
NFFT=2^nextpow2(L);
%transformada ponderada de Fourier para LPF
Y1=fft(LFP,NFFT)/L;
%obtenci?n de la mitad de los valores absolutos del espectro espejo de LFP
espectro1=2*abs(Y1(1:NFFT/2+1));
%transformada ponderada de Fourier para la envDemean
Y2=fft(envDemean,NFFT)/L;
%obtenci?n de la mitad de los valores absolutos del espectro espejo de la envDemean
espectro2=2*abs(Y2(1:NFFT/2+1));

%filtro de 100 puntos para suavizar ambos espectros
filter=ones(1,100)/100;
%convoluci?n de LFP con el filtro
espectro1=conv(espectro1,filter,'same');
%convoluci?n de la envDemean con el filtro
espectro2=conv(espectro2,filter,'same');

%obtenci?n del rango de frecuencias a trav?s de la frecuencia de sampleo
f=Fs/2*linspace(0,1,NFFT/2+1);
plot(f)


%% CALCULO DE LA POTENCIA DE LA envDemean
fr5=5;
tmp=abs(f-fr5);
[~, inx4]=min(tmp);
fr10=10;
tmp=abs(f-fr10);
[~, inx10]=min(tmp);
 
fr20=20;
tmp=abs(f-fr20);
[~, inx20]=min(tmp);
fr80=60;
tmp=abs(f-fr80);
[~, inx80]=min(tmp);
 
for n=inx4:inx10
    Area=(espectro2(n)+espectro2(n+1))/2*(f(n+1)-f(n));
end
GenvPow=sum(Area);
 
%% CALCULO DEL PICO MAXIMO DE FRECUENCIA
cerotes=zeros(length(espectro2(1:inx4-1)),1);%prepara un vector de 0 para aislar theta
newespectro2=vertcat(cerotes,espectro2(inx4:inx10));%concatena los ceros al rango de theta
[~,inxMaxE]=max(newespectro2);
GenvMax=f(inxMaxE);
 
cerotes=zeros(length(espectro1(1:inx20-1)),1);%prepara un vector de 0 para aislar gamma
newespectro1=vertcat(cerotes,espectro1(inx20:inx80));%concatena los ceros al rango de gamma
[~,inxMaxG]=max(newespectro1);
GMax=f(inxMaxG);
 
 
%% GRAFICAS
subplot(2,1,1)
plot(f,espectro1,f,espectro2)
axis([0 80 0 0.005]) %de 0 a 80 Hz y de 0 a 15e-3 mV2
subplot(2,1,2)
plot(Time,gamma,Time,envPeak,Time,theta)
axis([0 1 -0.35 0.35]) %de 0 a 1 s y de -0.1 a 0.1 mV
shg
 
%% Cross-correlation Gamma Envelope-Theta
%correlacion de las senyales
[acor,lag] = xcorr(envDemean,theta,'coeff');
 
 
[cor,I] = max((abs(acor)));%coeficiente de correlacion entre ambas senales 
Peaklag = lag(I);% Posicicion del pico maximo 
TimeDiff = Peaklag/Fs; % lag en s del pico maximo
Peaklagms = TimeDiff*1000; % Peak lag en ms.
 
plot(lag/Fs,acor)
axis([-0.5 0.5 -0.6 0.6])
shg
 
%% Autocorrelation of gamma signal and Rhythmicity
 
[acor,lag] = xcorr(gamma,500);      %correlacion de las senyales
Nacor = acor/max(acor); %Normaliza el maximo de autocorrelation
halfcor = Nacor(fix((length(acor)/2)):end); %mitad de la correlacion
halfTime = (lag(fix((length(lag)/2)):end))/Fs;%Time of the corresponding lag in sec
 
cor1 = halfcor; %redefine la variable halfcor
 
peaks = findpeaks(cor1);    %picos
troughs = findpeaks(-cor1); %valles
 
A = peaks(2) + 1;        %identificaci?n del primer pico
B = - troughs(1) + 1;  %identificaci?ndel segundo valle
 
Cr=(A-B)/(A+B);
 
%plot(halfTime,halfcor)
%axis([0 0.1 -1  1])
%shg
 
%% PAC package
Range = 1:1:100;
find_pac_shf(LFP, Fs, 'esc', LFP, Range, Range,'y',1,7,1000)
find_pac_shf(LFP, Fs, 'mi', LFP, Range, Range,'y',1,7,1000)
find_pac_shf(LFP, Fs, 'cfc', LFP, Range, Range,'y',1,7,1000)


