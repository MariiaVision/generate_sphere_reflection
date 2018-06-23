%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate scattering from a spherical 2 layer shell
%%%%%%%%%%%%%%% Mariia %%%%%% July 2017 %%%%%%%%


%%%%%%%%%%%%%% generate initial PULSE %%%%%%%%%%%%%%%
% pulse range:
f_pulse_0=160000; % [Hz] 
f_pulse_1=30000; % [Hz]
pulse_dur=0.0001; % pulse duration [sec]
fs=1000000; % [Hz]
NFFT=10001;

 pulse_generated=generate_chirp(f_pulse_0, f_pulse_1, pulse_dur, fs);
 pulse=pulse_generated';

s_analytic=hilbert(pulse);

% Frequencies of the pulse

s_freq=fft(s_analytic,NFFT);
s_ft_an=s_freq(1:NFFT/2+1);
freq_s=fs/2*linspace(0,1,NFFT/2+1);
f_step=freq_s(2)-freq_s(1); %based on the generated signal

figure();plot(freq_s, abs(s_ft_an)); title('original analytic signal'); xlabel('Frequency'), ylabel('Magnitude');

%%%%%%%%%%%%%%% Scattering calculation %%%%%%%%%%%%%%%%
%%%%%%% sphere's parameters %%%%%%%%%
    a=0.075; % radius [m]
    d=0.001; % thickness [m]
    r=1.5; % distance to the sphere [m]

% layer properties
p_1=1000; % kg/m^3 density of the outer liquid
p_2=2700;  % kg/m^3 density of the sphere

speed_long_1=1480; % m/sec longitudinal speed in the outer liquid (water)
speed_long_2= 6320; % m/sec longitudinal speed in the sphere material
speed_trans_2= 3130; % m/sec transversal speed in the sphere material


% 1st Lamé parameter
lame_2=p_2*(speed_long_2^2-2*speed_trans_2^2);
% shear modulus - 2nd lame parameter
mu_2=p_2*speed_trans_2^2;

p_3=1000; % density of the inner layer [kg/m^3]
speed_long_3=1480; % longitudinal speed in the inner layer  [m/sec]

% maximum number of modes
l=500;


%%%%%% calculate the ANALYTICAL REFLECTION
s_time=pulse;
[refl_time, refl_freq, time, freq, ff] = reflectionNumerical(s_time, s_freq, f_pulse_1, f_pulse_0, l, f_step, ...
    speed_long_1, speed_long_2, speed_long_3, speed_trans_2, r, a, d, p_1, p_2, p_3, NFFT, fs);


%%%% plot the results
figure();  plot(time*1000, refl_time); title('reflected signal'), xlabel('Time, msec'), ylabel('Amplitude');   %  the signal in time domain
figure(); plot( freq,fftshift(abs(refl_freq/max(refl_freq)))), title('Reflected signal'),xlabel('Frequency, Hz'), ylabel('Magnitude');   %the Fourier spectrum of x with zero frequency component


