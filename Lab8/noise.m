function [y_noisy] = noise(Eb_N0_dB,i,y,Nbps)
    for c = 1:size(y,1)
        Eb = var(y(c,:))/Nbps;
        EbN0 =  10^(Eb_N0_dB(i)/10); 
        N0 = Eb/EbN0;
        noise = sqrt(N0/2)*(randn(1,size(y,2))+1i*randn(1,size(y,2)));
        y_noisy(c,:) = y(c,:)+noise;
    end


end

