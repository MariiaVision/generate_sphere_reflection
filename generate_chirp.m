% The script is created to generate set of chirp-based pulses %
%%%             HWU, Ocean Systems Laboratory               %%%
%%%%%%%%%%%           Mariia Dmitrieva              %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pulse] =generate_chirp(f1, f2, dur, fs)

%dur = 0.001; % duration is one msec
%fs=1000000; % 1 MHz
t=0:1/fs:dur; % time sampling


%%%%%%%  WINDOWING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 1.7; % Gaussian
g = gausswin(length(t),alpha);

%Blackman
b=blackman(length(t));

a0=0.35875; 
a1=0.48829;  
a2=0.14128;
a3=0.01168;

b_manual=[];
N=length(t)-1;
%NT=3600;
for n=-N/2:N/2
b_manual=[b_manual,1-( a0 - a1*cos(2*pi*n/N) + a2*cos(2*pi*2*n/N))]; %-a3*cos(2*pi*n*3/N))];
end

%%% single chirp signals  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pulse_pre=chirp(t,f1,dur,f2, 'linear', -90); % -90 - set the right phase
pulse=pulse_pre.*b_manual;
pulse(1)=0;
pulse(end)=0;
                          % PLOT TIME & FREQUENCY %
%  time domain
figure();
plot(t, pulse,'b') 
xlabel('Time');


end