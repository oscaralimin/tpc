clc
clear all
close all

v1 = 10;        % Volt
duty = 0.75;    
len = 8;
T = 1e-6;       % sekunde
J = 1000;
samp_freq = J/T; % Hz

ind = 0.00001;  % H
cap = 0.000001; % F
res = 10;       % Ohm

% data sequence
input_raw = round(rand(1, len)); 

% create reference signal
ref_signal = zeros(1,J*len);
% duty cycle = percentage of uptime
for i = 1:len
    ref_signal(J*i-(J/2)-(duty*J/2):J*i-(J/2)+(duty*J/2)) = 1;
end

figure();
plot(ref_signal);

% For FSK
fsk_signal = [];
for i = 1:len
    fsk_signal = [fsk_signal generateFSK(duty, J, input_raw(i))];
end

figure();
plot(fsk_signal);

figure();
stem(input_raw);

% Calculate output
i_l = zeros(1, len*J);
v2 = zeros(1, len*J);
v2(1) = duty*v1;
i_l(1) = v2(1)/res;

tsl = 1/(samp_freq*ind);
tsc = 1/(samp_freq*cap);

for i = 2:J*len
    i_l(i) = i_l(i-1) - tsl*v2(i-1) + fsk_signal(i)*tsl*v1;
    v2(i) = v2(i-1) + tsc*i_l(i) - tsc/res*v2(i-1);
end

figure();
plot(v2);


function f = generateFSK(dutycycle, sfreq, bin)
   sfreq = fix(sfreq);
   sig = zeros(1,sfreq);
   if bin == 0
       sig((sfreq/2)-(dutycycle*sfreq/2):(sfreq/2)+(dutycycle*sfreq/2)) = 1;
   else
       sig((sfreq/4)-(dutycycle*sfreq/4):(sfreq/4)+(dutycycle*sfreq/4)) = 1;
       sig((sfreq*3/4)-(dutycycle*sfreq/4):(sfreq*3/4)+(dutycycle*sfreq/4)) = 1;
   end
   f = sig;
end






% duty = 0.25;
% ref_signal = zeros(1,samp_freq*length);
% for i = 1:length
%     ref_signal(samp_freq*i-(samp_freq/2)-(duty*samp_freq/2):samp_freq*i-(samp_freq/2)+1+(duty*samp_freq/2)) = 1;
% end
% plot(ref_signal)