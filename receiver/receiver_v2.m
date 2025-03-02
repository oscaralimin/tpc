% Function to simulate the receiver side. 
% 

function seq = receiver_v2(rx, matrix, sample_size, offset, s, len)
    % Downsample and perform MLSE on the received signal and compare to output matrix
    
    % Downsample received signal
    rx_temp = downsample(rx, sample_size/s, offset);
    
    err = zeros(1, 2^len);

    for i = 1:size(matrix, 1)
        % Downsample clean signal
        matrix_temp = downsample(matrix(i,:), sample_size/s, offset);
        
        % Calculate error
        err(i) = sum((rx_temp - matrix_temp).^2);
    end

    [~, index] = min(err);
    seq = int2bit(index-1, len)';
end

