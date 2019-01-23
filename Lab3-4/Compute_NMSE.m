function [sumNMSE] = Compute_NMSE(sumNMSE,c,channel_est_timefreq,channel_est,channels)
    if c == 3 
        NMSE = norm(channel_est_timefreq - fft(channels{1},size(channel_est_timefreq,2)),2)/norm(fft(channels{1},size(channel_est_timefreq,2)),2);
    elseif c == 4
        NMSE = norm(channel_est - fft(channels{2},size(channel_est,2)),2)/norm(fft(channels{2},size(channel_est,2)),2);
    else
        NMSE = 0;
    end
    sumNMSE = sumNMSE + NMSE;


end

