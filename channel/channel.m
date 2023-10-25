function rx = channel(sx)
    % add AWGN + modeling error (?)
    p = sum(sx.^2)/(length(sx));
    snr = 20;            % dB
    snr_lim = 10^(snr/10);
    const = sqrt(p)/sqrt(snr_lim);
    noise = const*randn(1, length(sx));
    rx = sx + noise;
end