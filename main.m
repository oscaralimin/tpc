clc
close all
clear all

%% Import other scripts
addpath("sender/","channel/","receiver/")

%% Parameters
flag_mod = 2;           % Flag for modulation scheme
                        % 0: Unmodulated
                        % 1: PWM
                        % 2: PSK
                        % 3: FSK

flag_lc = 0;            % Flag for line encoding
                        % 0: No line coding
                        % 1: Manchester coding

no_messages = 10000;      % # messages per SNR value

step_snr = 2;           % step size for SNR
min_snr = 0;          % min SNR value (dB)
max_snr = 20;           % max SNR value (dB)

downsamples = 4;
ds_offset = 100;

v1 = 10;        % Volt
duty = 0.75;    
len = 8;
T = 1e-6;      
sample_size = 1000;
samp_freq = sample_size/T;      % Hz
var = 0.2;      % PWM (duty1+duty2)/2 = duty       duty2 - duty1 = var

ind = 1e-5;     % H
cap = 1e-6;     % F
res = 15;       % Ohm

snr = min_snr:step_snr:max_snr;
ber = zeros(1,length(snr));




%% For loops
count = 1;
for i = min_snr:step_snr:max_snr        % SNR values in dB
    
    for k = 1:no_messages
        %% Sender
        [v2_apx, ~, ~, ~, send_seq] = sender(len, flag_lc, sample_size, flag_mod, duty, var, v1, cap, ind, res, samp_freq);
        sx = v2_apx;
        
        %% Channel
        [rx, sigma2] = channel(sx, duty*v1, i);
%         sum(sx.^2)/length(sx)
%         sum((sx-duty*v1).^2)/length(sx)
%         sum(rx.^2)/length(sx)
%         sigma2
%         pause()
        
        %% Receiver
        recv_seq = receiver(rx, len, flag_lc, sample_size, flag_mod, duty, var, v1, cap, ind, res, samp_freq, downsamples, ds_offset);
        
        %% Metric calculations
        ber(count) = ber(count) + biterr(send_seq, recv_seq);
%         a = biterr(send_seq, recv_seq, [], 'row-wise');
%         disp(send_seq)
%         disp(recv_seq)
%         disp(a);
%         
%         pause()
        
    end
    count = count + 1;
end

ber = ber/(no_messages*len);
figure
semilogy(min_snr:step_snr:max_snr,ber)
grid on
save workspace



