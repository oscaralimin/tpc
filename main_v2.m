clc
close all
clear all

%% Import other scripts
addpath("sender/","channel/","receiver/","buckConverter/", "ds_offset/")

%% Parameters
flag_mod = 2;           % Flag for modulation scheme
                        % 0: Unmodulated
                        % 1: PWM
                        % 2: PSK
                        % 3: FSK

flag_lc = 1;            % Flag for line encoding
                        % 0: No line coding
                        % 1: Manchester coding

no_messages = 50000;      % # messages per SNR value

step_snr = 1;           % step size for SNR
min_snr = 0;          % min SNR value (dB)
max_snr = 20;           % max SNR value (dB)

downsamples = 1;

v1 = 10;        % Volt
duty = 0.75;    
len = 4;
T = 1e-6;      
sample_size = 1000;
samp_freq = sample_size/T;      % Hz
var = 0.2;      % PWM (duty1+duty2)/2 = duty       duty2 - duty1 = var

ind = 1e-5;     % H
cap = 1e-6;     % F
res = 10;       % Ohm

snr = min_snr:step_snr:max_snr;
ber = zeros(1,length(snr));

%% Calculate every possible signal combinations
[v2_apx, ~, ~, ~] = buckConverter(flag_mod, flag_lc, duty, len, sample_size, ...
    samp_freq, v1, cap, ind, res, var);

%% Calculate optimal offset
if flag_lc == 0
    ds_offset = offset(v2_apx, len, downsamples);
elseif flag_lc == 1
    ds_offset = offset(v2_apx, len*2, downsamples);
end

count = 1;
%% For loop for snr and no_messages
for i = min_snr:step_snr:max_snr 
    for j = 1:no_messages
        %% Sender
        % Choose random sequence
        integer = randi([0, 2^len - 1]);
        send_seq = int2bit(integer, len)';
        sx = v2_apx(integer+1, :);
        
        %% Channel
        [rx, ~] = channel(sx, v1*duty, i);
    
        %% Receiver
        % Perform MLSE in the function
        recv_seq = receiver_v2(rx, v2_apx, sample_size, ds_offset, downsamples, len);
        
        %% Metric calculations
        ber(count) = ber(count) + biterr(send_seq, recv_seq);
    end
    count = count + 1;
end

ber = ber/(no_messages*len);
figure
semilogy(min_snr:step_snr:max_snr,ber)
grid on
save workspace







