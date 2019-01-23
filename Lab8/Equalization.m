function [equalized_signal] = Equalization(freq_signal,channel_est,Use_Multiple_groups_antenna,Use_different_distance,channels)
    if Use_Multiple_groups_antenna == true
        Channel_1_antenna = (abs(channel_est(1,:)).^2);
        Channel_2_antenna = 0;
        for j = 2:3
            Channel_2_antenna = Channel_2_antenna + (abs(channel_est(j,:)).^2);
        end
        Channel_3_antenna = 0;
        for j = 4:6
            Channel_3_antenna = Channel_3_antenna + (abs(channel_est(j,:)).^2);
        end
        equalized_signal = zeros(3,size(freq_signal,2),size(freq_signal,3));
        save_signal_1 = zeros(size(freq_signal,2),size(freq_signal,3));
        save_signal_2 = zeros(size(freq_signal,2),size(freq_signal,3));
        save_signal_3 = zeros(size(freq_signal,2),size(freq_signal,3));

        for j = 1:size(freq_signal, 2)
            save_signal_1(:,j) = ((conj(channel_est(1,:)).*freq_signal(1,:,j))./Channel_1_antenna).';
            save_signal_2(:,j) = (((conj(channel_est(2,:)).*freq_signal(2,:,j))./Channel_2_antenna).') + (((conj(channel_est(3,:)).*freq_signal(3,:,j))./Channel_2_antenna).');
            save_signal_3(:,j) = (((conj(channel_est(4,:)).*freq_signal(4,:,j))./Channel_3_antenna).') + (((conj(channel_est(5,:)).*freq_signal(5,:,j))./Channel_3_antenna).') + (((conj(channel_est(6,:)).*freq_signal(6,:,j))./Channel_3_antenna).');
        end
        equalized_signal(1,:,:) = save_signal_1;
        equalized_signal(2,:,:) = save_signal_2;
        equalized_signal(3,:,:) = save_signal_3;

    elseif Use_different_distance == true
        Channel_1_antenna = 0;
        Channel_2_antenna = 0;
        Channel_3_antenna = 0;
        for j = 1:3
            Channel_1_antenna = Channel_1_antenna + (abs(channel_est(j,:)).^2);
            Channel_2_antenna = Channel_2_antenna + (abs(channel_est(j+3,:)).^2);
            Channel_3_antenna = Channel_3_antenna + (abs(channel_est(j+6,:)).^2);
        end
        equalized_signal = zeros(3,size(freq_signal,2),size(freq_signal,3));
        save_signal_1 = zeros(size(freq_signal,2),size(freq_signal,3));
        save_signal_2 = zeros(size(freq_signal,2),size(freq_signal,3));
        save_signal_3 = zeros(size(freq_signal,2),size(freq_signal,3));

        for j = 1:size(freq_signal, 2)
            save_signal_1(:,j) = (((conj(channel_est(1,:)).*freq_signal(1,:,j))./Channel_1_antenna).') + (((conj(channel_est(2,:)).*freq_signal(2,:,j))./Channel_1_antenna).') + (((conj(channel_est(3,:)).*freq_signal(3,:,j))./Channel_1_antenna).');
            save_signal_2(:,j) = (((conj(channel_est(4,:)).*freq_signal(4,:,j))./Channel_2_antenna).') + (((conj(channel_est(5,:)).*freq_signal(5,:,j))./Channel_2_antenna).') + (((conj(channel_est(6,:)).*freq_signal(6,:,j))./Channel_2_antenna).');
            save_signal_3(:,j) = (((conj(channel_est(7,:)).*freq_signal(7,:,j))./Channel_3_antenna).') + (((conj(channel_est(8,:)).*freq_signal(8,:,j))./Channel_3_antenna).') + (((conj(channel_est(9,:)).*freq_signal(9,:,j))./Channel_3_antenna).');
        end
        equalized_signal(1,:,:) = save_signal_1;
        equalized_signal(2,:,:) = save_signal_2;
        equalized_signal(3,:,:) = save_signal_3;

    else
        Total_channel_est = 0;
        equalized_signal = zeros(1,size(freq_signal,2),size(freq_signal,3));
        save_signal = zeros(size(freq_signal,2),size(freq_signal,3));
        for j = 1:length(channels)
            Total_channel_est = Total_channel_est + (abs(channel_est(j,:)).^2);
        end
        for c = 1:length(channels)
           for j = 1:size(freq_signal, 2)
               save_signal(:,j) = ((conj(channel_est(c,:)).*freq_signal(c,:,j))./Total_channel_est).';
           end
           equalized_signal(1,:,:) = equalized_signal(1,:,:) + reshape(save_signal, [1 size(save_signal,1) size(save_signal,2)]);
        end
    end
end

