% Matlab function to simulate AWGN channel.
% Outputs the signal with added noise.

function [rx, sigma2] = channel(sx, snr, p)
    % add AWGN + modeling error (?)
    snr_lim = 10^(snr/10);
    const = sqrt(p)/sqrt(snr_lim);
    noise = const*randn(1, length(sx));
    
%     rx = awgn(sx-dc, snr, 'measured');
    rx = sx + noise;
    sigma2 = sum(noise.^2)/length(sx);
end