function seq = receiver(rx, len, flag_lc, sample_size, flag_mod, duty, var, v1, cap, ind, res, samp_freq, s, offset)
    lengthx = len;

    err = zeros(1, 2^lengthx);
    
    % Downsampling received signal
    rx_temp = downsample(rx, sample_size/s, offset);
%     rx_temp = rx;
%     figure()
%     clf
%     plot(rx)
%     hold on
%     plot([offset:sample_size/s:length(rx)],rx_temp,'o')

    signal_hypo = zeros(2^lengthx, sample_size*len);
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
        
        signal = [];
        if flag_mod == 0
            % For unmodulated
            for i = 1:len
                for j = 1:sample_size
                    signal = [signal input_raw(i)];
                end
            end
    
        elseif flag_mod == 1
            % For PWM
            for i = 1:len
                signal = [signal generatePWM(duty, sample_size, input_raw(i), var)];
            end
    
        elseif flag_mod == 2
            % For PSK
            for i = 1:len
                signal = [signal generatePSK(duty, sample_size,input_raw(i))];
            end
    
        elseif flag_mod == 3
            % For FSK
            for i = 1:len
                signal = [signal generateFSK(duty, sample_size, input_raw(i))];
            end
        end 
        
        [v2_apx, v2, i_l_apx, i_l] = buckConverter(signal, duty, len, sample_size, samp_freq, v1, cap, ind, res);
        signal_hypo(x+1,:) = v2_apx;
        % TODO
        % Downsampling
%         plot(v2_apx,'r')
        v2_apx = downsample(v2_apx, sample_size/s, offset);
%         plot([offset:sample_size/s:length(rx)],v2_apx,'rx')
%         pause()

        % Calculate error
        err(x+1) = sum((rx_temp - v2_apx).^2);
    end

%     save signal_hypotheses signal_hypo

    [~, index] = min(err);
    seq_temp = int2bit(index-1, lengthx)';

    if flag_lc == 1
        temp = [];
        for i = 1:length(seq_temp)
            if seq_temp(i) == 0
                temp = [temp 0 1];
            else
                temp = [temp 1 0];
            end
        end
        seq_temp = temp;
    end

    seq = seq_temp;
end

%% Modulation functions
function f = generateFSK(dutycycle, sfreq, bin)
   sfreq = fix(sfreq);
   sig = zeros(1,sfreq);
   if bin == 0
       sig((sfreq/2)-(dutycycle*sfreq/2):(sfreq/2)+(dutycycle*sfreq/2)) = 1;
   else
       sig(round((sfreq/4)-(dutycycle*sfreq/4)):round((sfreq/4)+(dutycycle*sfreq/4))) = 1;
       sig(round((sfreq*3/4)-(dutycycle*sfreq/4)):round((sfreq*3/4)+(dutycycle*sfreq/4))) = 1;
   end
   f = sig;
end

function f = generatePSK(dutycycle, sfreq, bin)
    sfreq = fix(sfreq);
    sig = zeros(1, sfreq);
    if bin == 0
        start_index = (sfreq/4) - (dutycycle*sfreq/2);
        end_index = (sfreq/4) + (dutycycle*sfreq/2);
        % if start index is less than index 1
        if start_index < 1
            temp = start_index;
            start_index = 1;
            sig(sfreq-abs(temp)+1:sfreq) = 1;
        end
        sig(start_index:end_index) = 1;
    else
        start_index = (sfreq*3/4) - (dutycycle*sfreq/2);
        end_index = (sfreq*3/4) + (dutycycle*sfreq/2);
        % if end_index is greater than sfreq
        if end_index > sfreq
            temp = end_index;
            end_index = sfreq;
            sig(1:temp-sfreq-1) = 1;
        end
        sig(start_index:end_index) = 1;
    end
    f = sig;
end

function f = generatePWM(dutycycle, sfreq, bin, var)
    sfreq = fix(sfreq);
    sig = zeros(1, sfreq);
    duty0 = dutycycle - var/2;
    duty1 = dutycycle + var/2;

    if bin == 0
        start_index = (sfreq/2) - (duty0*sfreq/2);
        end_index = (sfreq/2) + (duty0*sfreq/2);
        sig(start_index:end_index) = 1;
    else
        start_index = (sfreq/2) - (duty1*sfreq/2);
        end_index = (sfreq/2) + (duty1*sfreq/2);
        sig(start_index:end_index) = 1;
    end
    f = sig;
end

%% Buck converter
function [x, y, i, j] = buckConverter(signal, duty, len, sample_size, samp_freq, v1, cap, ind, res)
    % Output for approximation
    i_l_apx = zeros(1, len*sample_size);
    v2_apx = zeros(1, len*sample_size);
    v2_apx(1) = duty*v1;
    i_l_apx(1) = v2_apx(1)/res;
    
    tsl = 1/(samp_freq*ind);
    tsc = 1/(samp_freq*cap);
    
    % With approximation
    for i = 2:sample_size*len
        i_l_apx(i) = i_l_apx(i-1) - tsl*v2_apx(i-1) + signal(i)*tsl*v1;
        v2_apx(i) = v2_apx(i-1) + tsc*i_l_apx(i) - tsc/res*v2_apx(i-1);
    end
    
    % Output for w/o approximation
    i_l = zeros(1, len*sample_size);
    v2 = zeros(1, len*sample_size);
    v2(1) = duty*v1;
    i_l(1) = v2(1)/res;
    ts = 1/samp_freq;
    
    alpha = (cap*res*ts+ts^2) / (res*ts^2+cap*ind*res+ind*ts);
    beta = (cap*res*ts) / (res*ts^2+cap*ind*res+ind*ts);
    mu = (cap*res*ind+ind*ts) / (res*ts^2+cap*ind*res+ind*ts);
    gamma = (res*ts) / (cap*res+ts);
    delta = 1 - (ts/(cap*res+ts));
    
    % without approximation
    for i = 2:sample_size*len
        i_l(i) = mu*i_l(i-1) - beta*v2(i-1) + alpha*signal(i)*v1;
        v2(i) = gamma*i_l(i) + delta*v2(i-1);
    end
    
    x = v2_apx;
    y = v2;
    i = i_l_apx;
    j = i_l;    
end