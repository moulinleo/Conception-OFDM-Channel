function [y,channels] = Channel(H_new_LOS,H_new_NLOS,frame,Toa,channels,Time_acquisition,Use_different_distance,Use_Multiple_groups_antenna)

n=10;
if Use_different_distance == true
    channels = cell(9,1);
    % High spatial diversity
    channels{1} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{2} = reshape((H_new_LOS(6,6,6,:)),[1 n]);
    channels{3} = reshape((H_new_LOS(9,9,9,:)),[1 n]);
    % Middle spatial diversity
    channels{4} = reshape((H_new_LOS(9,9,9,:)),[1 n]);
    channels{5} = reshape((H_new_LOS(8,8,8,:)),[1 n]);
    channels{6} = reshape((H_new_LOS(10,10,10,:)),[1 n]);
    % Minimum spatial diversity
    channels{7} = reshape((H_new_LOS(4,4,4,:)),[1 n]);
    channels{8} = reshape((H_new_LOS(4,5,4,:)),[1 n]);
    channels{9} = reshape((H_new_LOS(4,4,5,:)),[1 n]);
elseif Use_Multiple_groups_antenna == true
    channels{1} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{2} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{3} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{4} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{5} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{6} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
else
    channels = cell(10,1);
    channels{1} = reshape((H_new_LOS(1,1,1,:)),[1 n]);
    channels{2} = reshape((H_new_LOS(2,2,2,:)),[1 n]);
    channels{3} = reshape((H_new_LOS(3,3,3,:)),[1 n]);
    channels{4} = reshape((H_new_LOS(4,4,4,:)),[1 n]);
    channels{5} = reshape((H_new_LOS(5,5,5,:)),[1 n]);
    channels{6} = reshape((H_new_LOS(6,6,6,:)),[1 n]);
    channels{7} = reshape((H_new_LOS(7,7,7,:)),[1 n]);
    channels{8} = reshape((H_new_LOS(8,8,8,:)),[1 n]);
    channels{9} = reshape((H_new_LOS(9,9,9,:)),[1 n]);
    channels{10} = reshape((H_new_LOS(10,10,10,:)),[1 n]);
end



for c = 1:length(channels)
    % Convolution with Channel
    y_conv = conv(channels{c},frame);
    y_conv = y_conv(1:end-(length(channels{c})-1));
    if Time_acquisition == true
        y_conv = [zeros(1,Toa) y_conv zeros(1,2000)];
    end 
    y(c,:) = y_conv;
end


end

