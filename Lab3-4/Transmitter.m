function [frame,preamble,preamble_time] = Transmitter(bit_tx,Nbps,modulation,N,cp_length)
symb = mapping(bit_tx',Nbps,modulation);

% S/P
stream = reshape(symb, [N length(symb)/N]);

% IFFT
stream_time = ifft(stream);

% CP
symb_ifft_cp = zeros(N + cp_length, length(symb)/N);
symb_ifft_cp(1:cp_length, :) = stream_time(end - cp_length + 1:end, :);
symb_ifft_cp(cp_length+1:end, :) = stream_time;

% P/S
sent = reshape(symb_ifft_cp,[1 numel(symb_ifft_cp)]);

% Preamble
A = 1;
preamble = datasample([-A A],N);
preamble_time = ifft(preamble);

pre_cp = zeros(1,2*cp_length+2*length(preamble_time));
pre_cp(1:cp_length) = preamble_time(end - cp_length + 1:end);
pre_cp(cp_length+1:2*cp_length) = preamble_time(end - cp_length + 1:end);
pre_cp(2*cp_length+1:end-length(preamble_time)) = preamble_time;
pre_cp(end-length(preamble_time)+1:end) = preamble_time;

frame = [pre_cp sent];

end

