% Matlab script for determining the minimum distance of buck converter
% signals as a function of subsampling (sample number and phase) and 
% sequence length

function f = offset(v2_mat, len, s)
    L = size(v2_mat, 2);
    J = L/len;
    s_offset = 0:J-1;

    min_dist_vec = zeros(1,length(s_offset));
    
    for jj = 1:length(s_offset)
    
        idx_s_cur = 1+s_offset(jj):J/s:L;
        signal_s = v2_mat(:,idx_s_cur);
            
        dist = [];
        for ii = 1:size(v2_mat,1)
            for kk = 1:size(v2_mat,1)
                if kk~=ii
                    dist = [dist sum(abs(signal_s(ii,:)-signal_s(kk,:)).^2)];
                end
            end
        end

        min_dist_vec(jj) = min(dist);
    end

    [~, idx] = max(min_dist_vec);
    f = idx - 1;
end
