function [y_noisy] = CFO_acquisition(y_noisy,Frequency_acquisition,Ts,delayEst,CFO_est,cp_length)
    if Frequency_acquisition == true
        for a =1:size(y_noisy,1)
            axis = 0:Ts:Ts*(length(y_noisy)-1);
            y_noisy(a,:) = y_noisy(a,:).*exp(-1j*(2*pi*CFO_est(delayEst+cp_length)*axis));
        end
    end
end

