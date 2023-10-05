clc
clear all
close all

v1 = 10;
duty = 0.25;
length = 8;
samp_freq = 1000;

% data sequence
input_raw = round(rand(1, length)); 

% create reference signal
ref_signal = zeros(1,samp_freq*length);
% duty cycle = percentage of uptime
for i = 1:length
    ref_signal(samp_freq*i-(samp_freq/2)-(duty*samp_freq/2):samp_freq*i-(samp_freq/2)+(duty*samp_freq/2)) = 1;
end

figure();
plot(ref_signal);

% For FSK
fsk_signal = [];
for i = 1:length
    fsk_signal = [fsk_signal generateFSK(duty, samp_freq, input_raw(i))];
end

figure();
plot(fsk_signal);

figure();
plot(input_raw);

function f = generateFSK(dutycycle, sfreq, bin)
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