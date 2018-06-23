%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numerical calculation of a Form Function based         %%%
%%% on the characteristics of a 2-layered spherical target %%%
%%% f_pulse_1 - lowest frequency to calculate from
%%% f_pulse_0 - highest frequency to calculate to
%%% l - number of modes
%%% f_step - frequency step
%%% speed_long_N - longitudinal speed in layers
%%% speed_trans_2 - transverse speed in the second layer
%%% p_N - density of the layers
%%%  lame_2, mu_2 - lame-parameters of the second layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mariia Dmitrieva June'17%

function [ff] = formfunction(f_pulse_1, f_pulse_0, l, f_step, speed_long_1, speed_long_2, speed_long_3, ...
    speed_trans_2, r_1, r_2, p_1, p_2, p_3, lame_2, mu_2)

%%% loops settings and vars

% frequency range:
f_start=f_pulse_1; %  kHz
f_end=f_pulse_0; % kHz
f_size=ceil((f_end-f_start)/f_step);

% mode range
l_start=0;
l_end=l;
l_size = l_end-l_start+1;

% set matrices for the loops
A_1=zeros(l_size, f_size);
I_matrix=zeros(l_size, f_size);
l_pos=0;

% mode loop
for l=l_start:l_end
    l_pos=l_pos+1;
    
    tic;
    f_pos=0;
    
    % frequency loop
    for f=f_start:f_step:f_end
        f_pos=f_pos+1;
        % wave numbers associated to the longitudinal and transversal sound speed
        k_L_1=2*pi*f/speed_long_1;
        k_L_2=2*pi*f/speed_long_2;
        k_L_3=2*pi*f/speed_long_3;
        k_T_2=2*pi*f/speed_trans_2;
        
        % longitudinal components
        x_L_1=k_L_1*r_1;
        x_L_2=k_L_2*r_2;
        
        y_L_1=k_L_2*r_1;
        y_L_2=k_L_3*r_2;
        
        %transversal components
        y_T_1=k_T_2*r_1;
        x_T_2=k_T_2*r_2;
        
        %%%%%%%%%%% find the cramer solution for the limitation equations %%%%%

        % form the matrix of coefficients
        
        A=  [ p_1/p_2*sbesselh(l, x_L_1), ...
             (lame_2*sbesselj(l,y_L_1) - 2*mu_2*ddbesselj(l, y_L_1))/(lame_2+2*mu_2), ...
             (lame_2*sbessely(l,y_L_1) - 2*mu_2*ddbessely(l, y_L_1))/(lame_2+2*mu_2), ...   
            -2*l*(l+1)/(y_T_1)^2 * (y_T_1*dbesselj(l, y_T_1)-sbesselj(l, y_T_1)) , ...
            -2*l*(l+1)/(y_T_1)^2 * (y_T_1*dbessely(l, y_T_1)-sbessely(l, y_T_1)) , ...
            0; ...
            
            x_L_1*dbesselh(l, x_L_1), ...
            y_L_1*dbesselj(l, y_L_1), ...
            y_L_1*dbessely(l, y_L_1), ...
            l*(l+1)*sbesselj(l, y_T_1)...
            l*(l+1)*sbessely(l, y_T_1)...
            0; ...
            
            0, ...
            2*(y_L_1*dbesselj(l, y_L_1)-sbesselj(l, y_L_1)), ...
            2*(y_L_1*dbessely(l, y_L_1)-sbessely(l, y_L_1)), ...
            y_T_1^2*ddbesselj(l, y_T_1)+(l+2)*(l-1)*sbesselj(l, y_T_1), ... % in Yan's thesis: (l+2)*(l+1)
            y_T_1^2*ddbessely(l, y_T_1)+(l+2)*(l-1)*sbessely(l, y_T_1), ...
            0; ...
            %%%%
            
            0, ...
            (lame_2*sbesselj(l,x_L_2) - 2*mu_2*ddbesselj(l, x_L_2))/(lame_2+2*mu_2), ...
            (lame_2*sbessely(l,x_L_2) - 2*mu_2*ddbessely(l, x_L_2))/(lame_2+2*mu_2), ...   
            -2*l*(l+1)/(x_T_2)^2 * (x_T_2*dbesselj(l, x_T_2)-sbesselj(l, x_T_2)) , ...
            -2*l*(l+1)/(x_T_2)^2 * (x_T_2*dbessely(l, x_T_2)-sbessely(l, x_T_2)) , ...
            p_3/p_2*sbesselj(l, y_L_2); ...
            
            0, ...
             x_L_2*dbesselj(l, x_L_2), ...
             x_L_2*dbessely(l, x_L_2), ...
             l*(l+1)*sbesselj(l, x_T_2)...
             l*(l+1)*sbessely(l, x_T_2)...
             y_L_2*dbesselj(l, y_L_2); ...
            
            0, ...
            2*(x_L_2*dbesselj(l, x_L_2)-sbesselj(l, x_L_2)), ...
            2*(x_L_2*dbessely(l, x_L_2)-sbessely(l, x_L_2)), ...
            x_T_2^2*ddbesselj(l, x_T_2)+(l+2)*(l-1)*sbesselj(l, x_T_2), ...  % in Yan's thesis: (l+2)*(l+1)
            x_T_2^2*ddbessely(l, x_T_2)+(l+2)*(l-1)*sbessely(l, x_T_2), ...
            0 ...
            ];
        
        % free components vector
        b=[(p_1/p_2)*sbesselj(l, x_L_1), x_L_1*dbesselj(l, x_L_1),0,0,0,0];
        
        [m, n] = size(b);
        
        %%% Soultion based on Cramer's rule
        
        % Cramer's Rule
        Ai = A;
        Ai(:,1) = b;
        A_det = (-1i^l*(2*l+1))*(det(Ai') / det(A')); % A matices are transformed to be on the right order
        
        A_1(l_pos,f_pos)=A_det;
        I_matrix(l_pos, f_pos)=1i^(l+3);
    end
        % break the loop when ff(30000)=NaN - limit for the mode number
    f_form_check=sum(A_1.*I_matrix); 
        if isnan(f_form_check(1))
        break
        end
end

%%%%%%% form function %%%%%%% 
ff=sum(A_1.*I_matrix); 

ff(isnan(ff))=0; % to eliminate NaN values which are coming up in case of high mode l for low frequencies beacuse of bessel functions

end
