function [sumBER] = Compute_BER(sumBER,equalized_signal,bit_tx,Nbps,modulation)
    for j = 1:size(equalized_signal,1)
        % P/S
        serial_signal = reshape(equalized_signal(j,:,:), [1, numel(equalized_signal(j,:,:))]);

        % Demapping
        bit_rx = demapping(serial_signal.',Nbps,modulation);

        % BER 
        diffs = abs(bit_tx' - bit_rx);
        errors = sum(diffs);
        BER = errors/length(bit_tx);
        sumBER(j) = sumBER(j) + BER;
    end
end

