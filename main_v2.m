% Script for simulating buck converter to be used for Talkative Power
% Conversion. At the end the Bit Error Rate will be calculated.
% Input parameters:
% - flag_mod        Flag for choosing modulation scheme
% - flag_lc         Flag to enable line encoding
%
% - no_messages     # of messages to be sent for each SNR value
% - step_snr        step size for SNR value
% - min_snr         minimum SNR value
% - max_snr         maximum SNR value
%
% - downsamples     # of samples for each symbol duration to be taken from 
%                   the output voltage
%
% - v1              Input voltage of the buck converter (V1)
% - duty            Duty cycle for operating the buck converter.
%                   V2 = duty cycle * V1
% - len             # of bits each message contains
% - T
% - sample_size     samples for each symbol
% - ind             inductance for the low-pass filter (Henry H)
% - cap             capacitance for the low-pass filter (Farad F)
% - res             resistance to act as load (Ohm)
%              
% - var             For PWM: difference between 2 duty cycle for bit 0 & 1
%
% Out:
% - ber             Bit error rate for each SNR value
% - wer             Word error rate

clc
close all
clear all

%% Import other scripts
addpath("sender/","channel/","receiver/","buckConverter/", "ds_offset/")

%% Parameters
flag_mod = 3;           % Flag for modulation scheme
                        % 0: Unmodulated
                        % 1: PWM
                        % 2: PSK
                        % 3: FSK

flag_lc = 0;            % Flag for line encoding
                        % 0: No line coding
                        % 1: Manchester coding

no_messages = 5000;    % # messages per SNR value

step_snr = 1;           % step size for SNR
min_snr = 0;            % min SNR value (dB)
max_snr = 20;           % max SNR value (dB)

downsamples = 8;        % # of samples to be taken for each symbol duration

v1 = 10;                % Input voltage (V1)
duty = 0.75;            % Duty cycle
len = 4;                % Bit length
T = 1e-6;      
sample_size = 1000;     % # samples for each symbol
ind = 1e-5;             % Henry
cap = 1e-6;             % Farrad
res = 10;               % Ohm

var = 0.2;              % PWM (duty1+duty2)/2 = duty       duty2 - duty1 = var

snr = min_snr:step_snr:max_snr;
ber = zeros(1,length(snr));
wer = zeros(1,length(snr));
samp_freq = sample_size/T;      % Hz

%% Calculate every possible signal combinations
[v2_apx, ~, ~, ~] = buckConverter(flag_mod, flag_lc, duty, len, sample_size, ...
    samp_freq, v1, cap, ind, res, var);
% Calculate noise variance for AWGN channel
p = sum(sum((v2_apx-v1*duty).^2))/(size(v2_apx,1)*size(v2_apx,2));

%% Calculate optimal offset
if flag_lc == 0
    [min_dis,ds_offset] = offset(v2_apx, len, downsamples);
elseif flag_lc == 1
    [min_dis,ds_offset] = offset(v2_apx, len*2, downsamples);
end

count = 1;
lb_snr = [];
snr_shift = 10*log10(downsamples/sample_size);
%% For loop for snr and no_messages
for i = min_snr:step_snr:max_snr 
    for j = 1:no_messages
        %% Sender
        % Choose random sequence
        integer = randi([0, 2^len - 1]);
        send_seq = int2bit(integer, len)';
        sx = v2_apx(integer+1, :);
        
        %% Channel
        [rx, ~] = channel(sx, i, p);
    
        %% Receiver
        % Perform MLSE in the function
        recv_seq = receiver_v2(rx, v2_apx, sample_size, ds_offset, downsamples, len);
        
        %% Metric calculations
        ber(count) = ber(count) + biterr(send_seq, recv_seq);
        if biterr(send_seq, recv_seq) > 0
            wer(count) = wer(count) + 1;
        end
    end
    lb_snr(count) = erfc( sqrt(min_dis*(10^( (i) /10) ) ) )/2;
    count = count + 1;
end

ber = ber/(no_messages*len);
wer = wer/(no_messages);
snr_shift = 10*log10(downsamples/sample_size);
figure
semilogy(snr,ber)
hold on
grid on
semilogy(snr+snr_shift,lb_snr, 'r')
axis([min_snr max_snr 1e-4 1])

figure
plot(snr, wer)
grid on
save workspace







