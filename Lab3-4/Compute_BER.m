function [sumBER] = Compute_BER(sumBER,bit_tx,bit_rx)
    diffs = abs(bit_tx' - bit_rx);
    errors = sum(diffs);
    BER = errors/length(bit_tx);
    sumBER = sumBER + BER;
end

