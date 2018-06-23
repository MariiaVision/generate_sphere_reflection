%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the function calculates reflection from a spherical
%%% 2-layered object with a given parameters
%%% output: refl_time - reflection in time domain
%%%         refl_freq - reflection in frequency domain
%%%         time - time axis
%%%         freq - freqiency axis
%%% input: 


%%%%%%%% Mariia %%%%%%%%%%%% July 2017 %%%%%%%%%%%%%%%%%

function [refl_time, refl_freq, time, freq, f_form] = reflectionNumerical(s, s_ft, f_pulse_1, f_pulse_0, l, f_step, ...
    speed_long_1, speed_long_2, speed_long_3, speed_trans_2, r, a, d, p_1, p_2, p_3,  NFFT, fs)

%%%%% calculate some parameters

% 1st Lamé parameter
lame_2=p_2*(speed_long_2^2-2*speed_trans_2^2);

% shear modulus - 2nd lame parameter
mu_2=p_2*speed_trans_2^2;

% radius of the layers' interface 
r_1=a; %+delta; %m
r_2=a-d; %m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        FORM FUNCTION          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_ft_an=s_ft(1:NFFT/2+1);


 ff=formfunction(f_pulse_1, f_pulse_0, l, f_step, speed_long_1, ...
     speed_long_2, speed_long_3, speed_trans_2, r_1, r_2, p_1, p_2, p_3, lame_2, mu_2);
 
% array for the form function
f_form=zeros(size(s_ft_an));

% calcualte where to insert the calculated values of f_form
n_start=(f_pulse_1/f_step);
n_end=(f_pulse_0/f_step);

% insert the calculated form function
f_form(n_start:n_end)=ff;
freq_s=fs/2*linspace(0,1,NFFT/2+1);

%%%%%% THE RESPONSE  %%%%%%%%%%%%%%%%%

k_L=2.*pi.*freq_s./speed_long_1;
k_L(k_L==0)=1;

% calculation of the response
[refl_time, refl_freq]=backscattering(k_L, r, s_ft_an, f_form, NFFT, fs, a);

freq=fs/2*linspace(-1,1,NFFT);
time=(0:length(refl_time)-1)*1/fs;
end