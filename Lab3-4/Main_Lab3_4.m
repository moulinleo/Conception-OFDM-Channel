%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      LAB 3-4     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    Multiple antenna receiver     %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Moulin Léo - Rouvroy Alexis - Schoone Rudy  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all

%% Load channel

data = load('Stat_channel');
h_NLOS_stat = data.h_NLOS_stat;
h_LOS_stat = data.h_LOS_stat;

%% Initialize parameters

Nbps = 4;
nb_bits = 8*4096;
bit_tx = randi([0,1],[1 nb_bits]);
Fs = 20e6;
modulation = 'qam';
N = 64;
cp_length = 16;
estimation = true;
Eb_N0_dB = -5:25;
channels = cell(4,1);
ber_vec = zeros(length(channels),length(Eb_N0_dB));
nmse_vec = zeros(length(channels),length(Eb_N0_dB));
essais = 50;

%% TRANSMITTER

[frame,preamble,preamble_time] = Transmitter(bit_tx,Nbps,modulation,N,cp_length);
%% Convolution with the New CHANNEL

channels{1} = h_LOS_stat;
channels{2} = h_NLOS_stat;
channels{3} = h_LOS_stat;
channels{4} = h_NLOS_stat;

for c = 1:length(channels)
    
    % Convolution with Channel
    y = conv(channels{c},frame);
    y= y(1:end-(length(channels{c})-1));
    
    for i = 1:length(Eb_N0_dB)
        sumBER = 0;
        sumNMSE = 0;
        for tries = 1:essais
            %% AWG Noise
            y_noisy = noise(Eb_N0_dB,i,y,Nbps);

            %% Channel estimation
            [y_noisy_nopre,channel_est_timefreq,channel_est] = Channel_estimation(y_noisy,cp_length,N,preamble,preamble_time,estimation,c,channels);

            %% Receiver
            bit_rx = Equalization(y_noisy_nopre,channel_est,channel_est_timefreq,cp_length,N,channels,c,Nbps,modulation);

            %% BER
            sumBER = Compute_BER(sumBER,bit_tx,bit_rx);
            
            %% NMSE
            sumNMSE = Compute_NMSE(sumNMSE,c,channel_est_timefreq,channel_est,channels);
        end    
        ber_vec(i,c) = sumBER/essais;
        nmse_vec(i,c) = sumNMSE/essais;  
    end
end

%% Plot BER

figure
semilogy(Eb_N0_dB,ber_vec(:,1),'r');
hold on;
semilogy(Eb_N0_dB,ber_vec(:,2),'b');
semilogy(Eb_N0_dB,ber_vec(:,3),'r--');
semilogy(Eb_N0_dB,ber_vec(:,4),'b--');

xlabel('Eb/N0 (dB)')
ylabel('BER')
grid on;
legend('NLOS','LOS','NLOS estimated channel','LOS estimated channel')
figure;
grid on;
plot(Eb_N0_dB,nmse_vec(:,3));
hold on;
plot(Eb_N0_dB,nmse_vec(:,4));
grid on;
xlabel('Eb/N0 (dB)')
ylabel('NMSE')
legend('Average in Time Domain','Average in Frequency Domain')

