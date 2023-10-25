function seq = receiver(rx)
    lengthx = length(rx);
    if flag_lc == 1
        lengthx = len/2;
    else
        lengthx = len;
    end

    err = zeros(1, 2^lengthx - 1);

    for x = 0:2^lengthx - 1
        input_raw = (int2bit(x,lengthx))';

        if flag_lc == 1
            temp = [];
            for i = 1:length(input_raw)
                if input_raw(i) == 0
                    temp = [temp 0 1];
                else
                    temp = [temp 1 0];
                end
            end
            input_raw = temp;
            len = length(input_raw);
        end
    end
end
