%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      LAB 8     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    Multiple antenna receiver     %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Moulin Léo - Rouvroy Alexis - Schoone Rudy  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all

%% Load new impulse responses
% 
H_new = load('Channels_new.mat');
H_new_LOS = H_new.H_new_LOS;
H_new_NLOS = H_new.H_new_NLOS;

%% Initialize parameters

Nbps = 4;
nb_bits = 4*4096;
bit_tx = randi([0,1],[1 nb_bits]);
Fs = 200e6;
Ts = 1/Fs;
modulation = 'qam';
N = 64;
cp_length = 16;
Time_acquisition = true;
Frequency_acquisition = false;
Use_Multiple_groups_antenna = true;
Use_different_distance = false;
CFO = 10;
phase_error = 0;
Eb_N0_dB = -5:25;
essais = 15;
channels = cell(6,1);
ber_vec = zeros(length(Eb_N0_dB));
nmse_vec = zeros(length(channels),length(Eb_N0_dB));
CFO_est =0;
delayEst=0;
%% TRANSMITTER ONE LINE

[frame,preamble,preamble_time] = Transmitter(bit_tx,Nbps,modulation,N,cp_length);
%% Add Time delay
if Time_acquisition == true
    Toa = randi(round(length(frame)*0.1));
else
    Toa = 0;
end


%% Convolution with the New CHANNEL
[y,channels] = Channel(H_new_LOS,H_new_NLOS,frame,Toa,channels,Time_acquisition,Use_different_distance,Use_Multiple_groups_antenna);

for i = 1:length(Eb_N0_dB)
    sumBER = zeros(3,1);
    sumNMSE = 0;
    for tries = 1:essais
        %% AWG Noise
        y_noisy = noise(Eb_N0_dB,i,y,Nbps);
        
        %% Add CFO
        y_noisy = Add_CFO(y_noisy,Frequency_acquisition,Ts,CFO,phase_error);

        %% RECEIVER

         %Estimate TOA
        [y_noisy,CFO_est,delayEst] = TOA_acquisition(y_noisy,frame,Time_acquisition,cp_length,CFO_est,delayEst);

        % Frequency acquisition
        y_noisy = CFO_acquisition(y_noisy,Frequency_acquisition,Ts,delayEst,CFO_est,cp_length);

        % Channel estimation + start of receiver line
        [freq_signal,channel_est] = Channel_estimation(y_noisy,cp_length,N,preamble,preamble_time);
        
        %Equalization        
        equalized_signal = Equalization(freq_signal,channel_est,Use_Multiple_groups_antenna,Use_different_distance,channels);
        
        % End of receiver line + BER computation
        sumBER = Compute_BER(sumBER,equalized_signal,bit_tx,Nbps,modulation);
    end
    for k=1:size(equalized_signal,1)
        ber_vec(i,k) = sumBER(k)/essais;
    end
end

%% Plot BER

if Use_Multiple_groups_antenna == true || Use_different_distance == true
    figure
    semilogy(Eb_N0_dB,ber_vec(:,1),'r');
    hold on;
    semilogy(Eb_N0_dB,ber_vec(:,2),'b');
    semilogy(Eb_N0_dB,ber_vec(:,3),'g');
    grid on
    xlabel('Eb/N0 (dB)')
    ylabel('BER')
    if Use_Multiple_groups_antenna == true
        legend('1 antenna at the receiver','2 antenna at the receiver','3 antenna at the receiver')
        title('BER performance as a function of the number of antennas')
    else
        legend('3 antenna at the receiver at pos (1,1,1);(6,6,6);(9,9,9)','3 antenna at the receiver at pos (8,8,8);(9,9,9);(10,10,10)','3 antenna at the receiver at pos (4,4,4);(4,5,4);(4,4,5)')
        title('Benefit of spatial diversity on BER')
    end
else
    figure
    semilogy(Eb_N0_dB,ber_vec(:,1),'r');
    grid on
    xlabel('Eb/N0 (dB)')
    ylabel('BER')
    legend([num2str(length(channels)) ' antenna at the receiver'])
end
