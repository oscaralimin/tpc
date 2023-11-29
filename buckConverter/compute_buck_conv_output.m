function [sig_v2_de] = compute_buck_conv_output(data,flag_mod,V1,R,C,L,T,delta,J)
%
% Matlab function for calulating the output voltage v2(t) of the buck 
% converter for different modulation schemes.
%
% In:
% - data       Data sequence vector (1xNp) containing transmit bits {0,1}
% - flag_mod   Flag for choosing the modulation scheme
%              0: unmodulated, 1: PWM, 2: PSK, 3: FSK
% - V1         Input voltage of buck converter in V
% - R          Resistor value in Ohm
% - C          Capacitor value in F
% - L          Inductor value in H
% - T          Bit duration in s
% - delta      Duty cycle (0 < delta <= 1)
% - J          Oversampling factor
%
% Out:
% - sig_v2_de  Output signal v2(t) calculated according via diff. equation
%
% J. Mietzner  
% Nov. 17, 2023
%
% *************************************************************************

% Preliminaries
Np = length(data);                     % Length of data bit sequence
dt = T/J;                              % Samping interval
t = [0:Np*J-1]*dt;                     % Time vector

% Pre-factors for solution of differential equation
a = 1/(2*C*R);
b_i = -sqrt( 1/(L*C) - a^2 );

% Initialize matrix for output signal (solution of diff. equation)
v2_indiv = zeros(Np+1,Np*J);  
% Contribution from steady state (solution of diff. equation)
v2_indiv(1,:) = delta*V1 * exp(-a*t) .* ( sin(b_i*t)/(C*R*b_i) - a*sin(b_i*t)/b_i+ cos(b_i*t) );

for ii = 2:size(v2_indiv,1)
    
    % Initializations
    v2_pulse_cur_alpha_zp = zeros(1,size(v2_indiv,2));
    v2_pulse_cur_beta_zp = zeros(1,size(v2_indiv,2));
    start_idx_nT = (ii-2)*J + 1;
    length_sig = size(v2_indiv,2) - start_idx_nT + 1;
    
    % Signal part shifted by alpha_n
    if flag_mod == 0
        % Pulse start for unmodulated case
        alpha_n = T/2 - (delta*T)/2;
    elseif flag_mod == 1
        % Pulse start for PWM case
        delta_short = delta - 0.1;
        delta_long = delta + 0.1;
        if data(ii-1) == 0
            delta_act = delta_short;
        else
            delta_act = delta_long;
        end
        alpha_n = T/2 - (delta_act*T)/2;
    elseif flag_mod == 2
        % Pulse start for PSK case
        if delta < 0.5
            if data(ii-1) == 0
                alpha_n = T/4 - (delta*T)/2;
            else
                alpha_n = 3*T/4 - (delta*T)/2;
            end
        else
            if data(ii-1) == 0
                alpha_n = 3*T/4 + (delta*T)/2;
            else
                alpha_n = T/4 + (delta*T)/2;
            end
        end
    elseif flag_mod == 3
        % Pulse start(s) for FSK case
        if data(ii-1) == 0 % - Only one FSK pulse
            alpha_n = T/2 - (delta*T)/2;
        else % - Two FSK pulses
            alpha_n = T/4 - (delta*T)/4;
            alpha_n2 = 3*T/4 - (delta*T)/4;
        end        
    end
    t_s_ov = max(t-alpha_n,0);
    v2_pulse_cur_alpha = V1 * ( 1 - exp(-a*t_s_ov).*cos(b_i*t_s_ov)  - a/b_i*exp(-a*t_s_ov).*sin(b_i*t_s_ov) );
    v2_pulse_cur_alpha_zp(start_idx_nT:end) = v2_pulse_cur_alpha(1:length_sig);
    
    % Account for pulse start at t=0 for PSK with delta >= 0.5
    if 0 % - currently disabled 
    if (flag_mod == 2)&&(delta >= 0.5)&&(ii == 2)
        v2_pulse_cur_alpha = V1 * ( 1 - exp(-a*t).*cos(b_i*t)  - a/b_i*exp(-a*t).*sin(b_i*t) );
        v2_pulse_cur_alpha_zp(start_idx_nT:end) = v2_pulse_cur_alpha_zp(start_idx_nT:end) + v2_pulse_cur_alpha(1:length_sig);
    end
    end    
    
    % Account for second pulse in FSK scheme for data bit = 1
    if (flag_mod == 3)&&(data(ii-1) == 1)
        t_s_ov = max(t-alpha_n2,0);
        v2_pulse_cur_alpha = V1 * ( 1 - exp(-a*t_s_ov).*cos(b_i*t_s_ov)  - a/b_i*exp(-a*t_s_ov).*sin(b_i*t_s_ov) );
        v2_pulse_cur_alpha_zp(start_idx_nT:end) = v2_pulse_cur_alpha_zp(start_idx_nT:end) + v2_pulse_cur_alpha(1:length_sig);
    end    
 
    % Signal part shifted by beta_n
    if flag_mod == 0
        % Pulse end for unmodulated case
        beta_n = T/2 + (delta*T)/2;
    elseif flag_mod == 1
        % Pulse end for PWM case
        beta_n = T/2 + (delta_act*T)/2;
    elseif flag_mod == 2
        % Pulse end for PSK case
        if delta < 0.5
            if data(ii-1) == 0
                beta_n = T/4 + (delta*T)/2;
            else
                beta_n = 3*T/4 + (delta*T)/2;
            end
        else
            if data(ii-1) == 0
                beta_n = 3*T/4 - (delta*T)/2;
            else
                beta_n = T/4 - (delta*T)/2;
            end
        end
    elseif flag_mod == 3
        % Pulse ends(s) for FSK case
        if data(ii-1) == 0 % - Only one FSK pulse
            beta_n = T/2 + (delta*T)/2;
        else % - Two FSK pulses
            beta_n = T/4 + (delta*T)/4;
            beta_n2 = 3*T/4 + (delta*T)/4;
        end
    end
    t_s_ov = max(t-beta_n,0);
    v2_pulse_cur_beta = -V1 * ( 1 - exp(-a*t_s_ov).*cos(b_i*t_s_ov)  - a/b_i*exp(-a*t_s_ov).*sin(b_i*t_s_ov) );
    v2_pulse_cur_beta_zp(start_idx_nT:end) = v2_pulse_cur_beta(1:length_sig);  
    
    % Account for second pulse in FSK scheme for data bit = 1
    if (flag_mod == 3)&&(data(ii-1) == 1)
        t_s_ov = max(t-beta_n2,0);
        v2_pulse_cur_beta = -V1 * ( 1 - exp(-a*t_s_ov).*cos(b_i*t_s_ov)  - a/b_i*exp(-a*t_s_ov).*sin(b_i*t_s_ov) );
        v2_pulse_cur_beta_zp(start_idx_nT:end) = v2_pulse_cur_beta_zp(start_idx_nT:end) + v2_pulse_cur_beta(1:length_sig);  
    end    
    
    % Overall contribution of pulse number (ii-1)
    v2_indiv(ii,:) = v2_pulse_cur_alpha_zp + v2_pulse_cur_beta_zp;
    
end

% Final output voltage
sig_v2_de = sum(v2_indiv,1);
