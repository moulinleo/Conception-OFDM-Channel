function [y_noisy] = noise(Eb_N0_dB,i,y,Nbps)
    % Add noise to signal
    Eb = var(y)/Nbps;
    EbN0 =  10^(Eb_N0_dB(i)/10); 
    N0 = Eb/EbN0;
    noise = sqrt(N0/2)*(randn(1,length(y))+1i*randn(1,length(y)));
    y_noisy = y + noise;
end

