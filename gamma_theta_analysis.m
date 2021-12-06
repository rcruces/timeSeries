clear
clc

%% we read data
data = readtable('basal.csv');
%sampling rate
srate = 1/(data.times(2)-data.times(1));

%% band pass filter
%we apply a band pass filter to keep only those frequencies between 1 and
%120 Hz
[signal_filt,filt]=bandpass(data.signal,[3 120],srate);

figure(1), clf

plot(data.times,data.signal) 
hold on 
plot(data.times,signal_filt)
hold off
xlabel('Time (s)')
ylabel('mV')
title('Recorded signal')
legend('Unfiltered signal','Filtered signal')

%uncomment for plot
%bandpass(data.signal,[1 120],srate)

%% frequency anañysis
fourierCoefs = fft(signal_filt) / length(data.times);
hz = linspace(0,srate/2,floor(length(data.times)/2)+1);

figure(2), clf
plot(hz, abs(fourierCoefs(1:length(data.times)/2+1)));
title('Frequency domain')
xlabel('Frequencies (Hz)')
ylabel('Power')
xlim([0,120])

%% Waves extraction
%Theta wave
theta = hz>3 & hz<8;
theta_coefs = fourierCoefs(theta);
theta_power = mean(abs(theta_coefs)) %Mean theta power
hz_theta = linspace(3, 8, length(theta_coefs));
max_theta = max(abs(theta_coefs)); %Maximum theta power
theta_index = find(abs(theta_coefs)==max_theta);
theta_peak = hz_theta(theta_index) %Gamma peak frequency
total_theta_power = sum(abs(theta_coefs))

figure(3), clf
plot(hz_theta, abs(theta_coefs)*2)
title('Frequency domain (theta wave)')
xlabel('Frequencies (Hz)')
ylabel('Amplitude')
xlim([3,8])


%Gamma wave
gamma = hz>20 & hz<80;
gamma_coefs = fourierCoefs(gamma);
hz_gamma = linspace(20, 40, length(gamma_coefs));
gamma_power_mean = mean(abs(gamma_coefs)) %Mean gamma power
max_gamma = max(abs(gamma_coefs)); %Maximum gamma power
gamma_index = find(abs(gamma_coefs)==max_gamma);
gamma_peak = hz_gamma(gamma_index) %Gamma peak frequency
total_gamma_power = sum(abs(gamma_coefs))

figure(4), clf
plot(hz_gamma, abs(gamma_coefs)*2)
title('Frequency domain (gamma wave)')
xlabel('Frequencies (Hz)')
ylabel('Amplitude')
xlim([20,40])

%ratio of gamma and theta amplitudes
ratio = mean(abs(gamma_coefs))/mean(abs(theta_coefs)) 

%plot of gamma, theta and signal
[theta_signal,theta_filt]=bandpass(data.signal,[3 8],srate);
[gamma_signal,gamma_filt]=bandpass(data.signal,[40 120],srate);

figure(5), clf

plot(data.times, signal_filt) 
hold on 
plot(data.times, theta_signal) 
plot(data.times, gamma_signal)
hold off
xlabel('Time (s)')
ylabel('mV')
title('Signals')
legend('Signal','Theta wave','Gamma wave')

%% Envelope signal of gamma

[gamma_upper,gamma_lower] = envelope(gamma_signal, 1000,'peak');

fourierCoefs_envelope = fft(gamma_upper) / length(data.times);

figure(6), clf
plot(hz, abs(fourierCoefs_envelope(1:length(data.times)/2+1)));
title('Frequency domain (envelope signal)')
xlabel('Frequencies (Hz)')
ylabel('Power')
xlim([0,3])

figure(7), clf

plot(data.times, gamma_signal)
hold on 
plot(data.times, gamma_upper)
hold off
xlabel('Time (s)')
ylabel('mV')
title('Signals')
legend('Gamma wave','Envelope wave')

%% Autocorrelation analysis
autocorr(gamma_upper)

