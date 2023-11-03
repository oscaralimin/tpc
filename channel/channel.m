function [rx, sigma2] = channel(sx, dc, snr)
    % add AWGN + modeling error (?)
    p = sum((sx-dc).^2)/(length(sx));
    snr_lim = 10^(snr/10);
    const = sqrt(p)/sqrt(snr_lim);
    noise = const*randn(1, length(sx));
    
%     rx = awgn(sx-dc, snr, 'measured');
    rx = sx + noise;
    sigma2 = sum(noise.^2)/length(sx);
end