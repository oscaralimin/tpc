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

flag_lc = 1;            % Flag for line encoding
                        % 0: No line coding
                        % 1: Manchester coding

no_messages = 10;     % # messages per SNR value                    

v1 = 10;        % Volt
duty = 0.75;    
len = 8;
T = 1e-6;      
sample_size = 1000;
samp_freq = sample_size/T;      % Hz
var = 0.2;      % PWM (duty1+duty2)/2 = duty       duty2 - duty1 = var

ind = 1e-5;     % H
cap = 1e-6;     % F
res = 10;       % Ohm

%% For loops
for i = 0:0.2:20        % SNR values in dB
    for k = 1:no_messages
%% Sender
[v2_apx, v2, i_l_apx, i_l] = sender(len, flag_lc, sample_size, flag_mod, duty, var, v1, cap, ind, res, samp_freq);
sx = v2_apx;

%% Channel
rx = channel(sx);

%% Receiver


%% Metric calculations

    end
end