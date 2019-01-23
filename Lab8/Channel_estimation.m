function [freq_signal,channel_est] = Channel_estimation(y_noisy,cp_length,N,preamble,preamble_time)
    for c = 1:size(y_noisy,1)
        % Preamble & CP Removal
        preamble_out1 = y_noisy(c,2*cp_length+1:2*(cp_length)+N);
        preamble_out2 = y_noisy(c,2*(cp_length)+N+1:2*(cp_length+N));
        y_noisy_save= y_noisy(c,2*(cp_length+N)+1:end);  
        y_noisy_nopre(c,:) = y_noisy_save;


        channel_est_time = preamble_time\((preamble_out1 + preamble_out2)/2);
        channel_est_timefreq = fft(channel_est_time(1,:));

        % FFT Preamble
        preamble_out_freq1 = fft(preamble_out1);
        preamble_out_freq2 = fft(preamble_out2);

        channel_est1 = preamble_out_freq1./preamble;
        channel_est2 = preamble_out_freq2./preamble;
        channel_est(c,:) = (channel_est1 + channel_est2)/2;

        % S/P
        y_rec(c,:,:) = reshape(y_noisy_nopre(c,:), [1, N + cp_length, length(y_noisy_nopre(c,:))/(N + cp_length)]);

        % CP Removal
        signal_no_cp(c,:,:) = y_rec(c,cp_length + 1:end, :);

        % FFT
        freq_signal(c,:,:) = fft(signal_no_cp(c,:,:));
    end
end

