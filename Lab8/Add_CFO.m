function [y_noisy] = Add_CFO(y_noisy,Frequency_acquisition,Ts,CFO,phase_error)
    if Frequency_acquisition == true
        for a=1:size(y_noisy,1)
            axis = 0:Ts:Ts*(length(y_noisy)-1);
            y_noisy(a,:) = y_noisy(a,:).*exp(1j*(2*pi*CFO*axis+phase_error));
        end
    end
end

