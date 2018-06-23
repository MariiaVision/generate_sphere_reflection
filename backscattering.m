%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the function calculates backscattered signal based on 
%%%                                       the form function
%%%
%%% output: refl - reflection in time domain
%%%         refl_ft - reflection in frequency domain
%%% input:  k_L - arrays with wave numbers for the set of frequencies
%%%         r - distance from the target center to the pulse source 
%%%         s_ft_an - orinial pulse in FD 
%%%         f_form - form function 
%%%         NFFT - number of samples in FD 
%%%         fs - sampling frequency 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mariia Dmitrieva - June'17 %%%


function [refl, refl_ft]=backscattering(k_L, r, s_ft_an, f_form, NFFT, fs, a)

% response in frequency domain for the analytic signal
res_ft=(exp(-2i.*k_L'.*r))./(abs(k_L').*r.^2).*s_ft_an.*f_form;

res_ft(1)=0; % to be sure that it is zero


% echo in time domain with or without circle shift on NFFT/2;
refl=circshift(real(ifft(res_ft, NFFT)), floor(NFFT/2));
refl=real(ifft(res_ft, NFFT));

refl_ft=fft(refl,NFFT);

end