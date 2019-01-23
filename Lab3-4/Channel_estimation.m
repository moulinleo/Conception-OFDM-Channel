function [y_noisy_nopre,channel_est_timefreq,channel_est] = Channel_estimation(y_noisy,cp_length,N,preamble,preamble_time,estimation,c,channels)
        % Preamble & CP Removal
        preamble_out1 = y_noisy(2*cp_length+1:2*(cp_length)+N);
        preamble_out2 = y_noisy(2*(cp_length)+N+1:2*(cp_length+N));
        y_noisy_nopre = y_noisy(2*(cp_length+N)+1:end);  

        channel_est_time = preamble_time\((preamble_out1 + preamble_out2)/2);
        channel_est_timefreq = fft(channel_est_time(1,:));

        % FFT Preamble
        preamble_out_freq1 = fft(preamble_out1);
        preamble_out_freq2 = fft(preamble_out2);

        channel_est1 = preamble_out_freq1./preamble;
        channel_est2 = preamble_out_freq2./preamble;
        channel_est = (channel_est1 + channel_est2)/2;
end

