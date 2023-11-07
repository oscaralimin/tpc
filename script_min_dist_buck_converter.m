% Matlab script for determining the minimum distance of buck converter
% signals as a function of subsampling (sample number and phase) and 
% sequence length
% 
close all

% Plot option
offset_plot = 0;              % Choose sampling offset for example plots

% Load example signals 
%load signal_hypotheses_L4     % PSK modulation without line coding, 4 bits
%n_bits = 4;
load signal_hypotheses_L8     % PSK modulation without line coding, 8 bits
n_bits = 8;

% Sampling process and assessment of minimum distance
L = size(signal_hypo,2);
J = L/n_bits;
S = 2;
s_offset = [0:J-1];

min_dist_vec = zeros(1,length(s_offset));

for jj = 1:length(s_offset)

    idx_s_cur = [1+s_offset(jj):J/S:L];
    signal_s = signal_hypo(:,idx_s_cur);

%     % Test plot
%     if (jj-1)==offset_plot
%         figure
%         plot(signal_hypo.')
%         hold on
%         plot(idx_s_cur,signal_s.','ro')
%         grid on
%         title(['Length of bit sequence = ' num2str(n_bits) ', S = ' num2str(S) ', J = ' num2str(J)])
%     end
        
    dist = [];
    for ii = 1:size(signal_hypo,1)
        for kk = 1:size(signal_hypo,1)
            if kk~=ii
                dist = [dist sum(abs(signal_s(ii,:)-signal_s(kk,:)).^2)];
            end
        end
    end
    min_dist_vec(jj) = min(dist);
    
end

figure
plot(s_offset,min_dist_vec,'b-')
grid on
xlabel('Sample offset (# samples)')
ylabel('Minimum distance')
title(['Length of bit sequence = ' num2str(n_bits) ', S = ' num2str(S) ', J = ' num2str(J)])




