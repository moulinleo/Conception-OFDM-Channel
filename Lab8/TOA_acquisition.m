function [y_noisy,CFO_est,delayEst] = TOA_acquisition(y_noisy,frame,Time_acquisition,cp_length,CFO_est,delayEst)
     if Time_acquisition == true
        lengthCorr = round(length(frame)*0.1);
        correlation_vector = zeros(1,lengthCorr);
        somme = zeros(1,lengthCorr);

        for n = 1:lengthCorr
            for l = 0:64-1
                correlation_vector(n) =  correlation_vector(n) + y_noisy(1,n+l+64).* conj(y_noisy(1,n+l));
            end
            for l = 0: 2*64 -1
                somme(n) = somme (n) + abs(y_noisy(1,n+l)).^2;
            end
        end

        CFO_est = (phase(correlation_vector));

        correlation = abs(correlation_vector)./(somme);
        M = zeros(1,length(correlation));
        AvgSize = 32;

        for k = 1:length(correlation)-AvgSize
            M(k) = sum(abs(correlation(k:k+AvgSize)))/(AvgSize+1);
        end

        [ ~, delayEst] = max (M);

        if CFO_est(delayEst+cp_length)<0
            CFO_est = CFO_est + 2*pi;
        end
        if CFO_est(delayEst+cp_length)>3.14
            CFO_est = CFO_est - 2*pi;
        end

        CFO_est = CFO_est*1000;

        toa_est = delayEst-15;

        if(toa_est<=0)
           toa_est = 1; 
        end

        y_noisy = y_noisy(:,toa_est : toa_est + length(frame)-1);     
     end
end

