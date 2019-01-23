function [bit_rx] = Equalization(y_noisy_nopre,channel_est,channel_est_timefreq,cp_length,N,channels,c,Nbps,modulation)
        % S/P
        y_rec = reshape(y_noisy_nopre, [N + cp_length, length(y_noisy_nopre)/(N + cp_length)]);

        % CP Removal
        signal_no_cp = y_rec(cp_length + 1:end, :);

        % FFT
        freq_signal = fft(signal_no_cp);   
        
        % Equalization
        equalized_signal = zeros(size(freq_signal));
        if c == 1 || c == 2
            for j = 1:size(freq_signal, 2)
                equalized_signal(:,j) = freq_signal(:,j) ./ fft(channels{c}, size(freq_signal, 1)).';
            end
        else
            for j = 1:size(freq_signal, 2)
                equalized_signal(:,j) = freq_signal(:,j) ./ channel_est.';
            end
        end
        
        % P/S
        serial_signal = reshape(equalized_signal, [1, numel(equalized_signal)]);

        % Demapping
        bit_rx = demapping(serial_signal.',Nbps,modulation);
end

