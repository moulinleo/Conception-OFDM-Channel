%% LABO 5
clear all;
close all;
clc;

% Load impulse responses
data = load('h_wide_NLOS_stat.mat');
h_wide_NLOS_stat = data.h_wide_NLOS_stat;
data2 = load('h_wide_LOS_stat.mat');
h_wide_LOS_stat = data2.h_wide_LOS_stat;


%% PARAMETERS
Nbps = 4;
nb_bits = 8*4096;
bit_tx = randi([0,1],[1 nb_bits]);
Fs = 20e6;
modulation = 'qam';
N = 64;
cp_length = 16;
estimation = false;

%% TRANSMITTER
symb = mapping(bit_tx',Nbps,modulation);
% scatterplot(symb)

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
preamble = datasample([-1 1],N);
preamble_time = ifft(preamble);
power_signal = var(sent)/Nbps;
power_preamble = var(preamble_time)/Nbps;

% Same power preamble & signal
preamble_time = preamble_time*power_signal/power_preamble;

pre_cp = zeros(1,2*cp_length+2*length(preamble_time));
pre_cp(1:cp_length) = preamble_time(end - cp_length + 1:end);
pre_cp(cp_length+1:2*cp_length) = preamble_time(end - cp_length + 1:end);
pre_cp(2*cp_length+1:end-length(preamble_time)) = preamble_time;
pre_cp(end-length(preamble_time)+1:end) = preamble_time;

frame = [pre_cp sent];

% Time delay
Toa = randi(round(length(frame)*0.1));
Toa = 500;


%% CHANNEL
channels = cell(4,1);
channels{1} = h_wide_NLOS_stat;
channels{2} = h_wide_LOS_stat;
channels{3} = h_wide_NLOS_stat;
channels{4} = h_wide_LOS_stat;
Eb_N0_dB = -5:25;
ber_vec = zeros(length(channels),length(Eb_N0_dB));
nmse_vec = zeros(length(channels),length(Eb_N0_dB));
rmse_vec = zeros(length(channels),length(Eb_N0_dB));
essais = 1;

for c = 1:length(channels)
    % Convolution with Channel
    y = conv(channels{c},frame);
    y = y(1:end-5);
    y = [zeros(1,Toa) y zeros(1,2000)];
    for i = 1:length(Eb_N0_dB)
        sumBER = 0;
        sumNMSE = 0;
        sumRMSE = 0;
        for tries = 1:essais
            %% AWG Noise
            Eb = var(y)/Nbps;
            EbN0 =  10^(Eb_N0_dB(i)/10); 
            N0 = Eb/EbN0;
            noise = sqrt(N0/2)*(randn(1,length(y))+1i*randn(1,length(y)));
            y_noisy = y + noise ;
            

            %% RECEIVER
             % Estimate TOA
             if c == 1 || c == 2
                correlation_vector = zeros(1,round(length(frame)*0.1));
                
                %% Cross-correlation with preamble
                for ii = 1:round(length(frame)*0.1)
                    corr1 = [abs(preamble_time) abs(preamble_time)];
                    corr2 = abs(y_noisy(ii+1:ii+2*N));
                    correlation_vector (ii) = corr(corr1.',corr2.');
                end
                %% Autocorrelation between successive blocks
%                 for ii = 1:round(length(frame)*0.1)
%                     corr1 = abs(y_noisy(1+ii:N+ii));
%                     corr2 = abs(y_noisy(N+1+ii:2*N+ii));
%                     correlation_vector (ii) = corr(corr1.',corr2.');
%                 end
                [a in] = max(correlation_vector);
                dd = Toa + cp_length;
%                 if i == 1
%                     figure;
%                     plot(correlation_vector); hold on;
%                     line([dd dd],[min(correlation_vector) max(correlation_vector)],'Color','green','LineWidth',1);
%                     line([in in],[min(correlation_vector) max(correlation_vector)],'Color','red','LineWidth',1);
%                     title(['CrossCorrelation with preamble. EbN0 = ' num2str(Eb_N0_dB(i)) ' dB']);
%                     legend('CrossCorr','True TOA','Estimated TOA = max(CrossCorr)');
%                     xlabel('Time steps'); ylabel('Autocorrelation');
%                     grid on;
%                 end
                
                %Moving average window
                size_plateau = sum(correlation_vector>(max(correlation_vector)-0.1));
                peak = zeros(round(length(frame)*0.1)-size_plateau,1);
                for kk = 1:round(length(frame)*0.1)-size_plateau
                    peak(kk) = mean(correlation_vector(kk:kk+size_plateau));
                end
                [~,inx] = max(peak);
                %toa_est = (inx + round(0.5*size_plateau))-2*cp_length
                %toa_est = inx + round(0.5*size_plateau) + 5;
                toa_est = inx;
                %toa_Est = Toa;
                dd = toa_est;
                if i == 1
                    figure();
                    plot(peak); hold on;
                    %line([dd dd],[min(peak) max(peak)],'Color','green','LineWidth',1);
                    %line([toa_est toa_est],[min(peak) max(peak)],'Color','red','LineWidth',1);
                    plot(toa_est,max(peak), 'x');
                    title(['CrossCorrelation with preamble, EbN0 = ' num2str(Eb_N0_dB(i)) ' dB']);
                    %legend('CrossCorr','True TOA','Estimated TOA = max(CrossCorr)');
                    legend('CrossCorr','Estimated TOA');
                    xlabel('Time steps'); ylabel('CrossCorrelation');
                    grid on;
                end
                
                y_noisy = y_noisy(toa_est : toa_est + length(frame)-1);     
            
                
            else
                y_noisy = y_noisy(Toa+1 : Toa + length(frame)); 
            end
            
            % Preamble & CP Removal
            preamble_out1 = y_noisy(2*cp_length+1:2*(cp_length)+N);
            preamble_out2 = y_noisy(2*(cp_length)+N+1:2*(cp_length+N));
            y_noisy_nopre = y_noisy(2*(cp_length+N)+1:end);  
            
            channel_est_time = preamble_time\((preamble_out1 + preamble_out2)/2);
            channel_est_timefreq = fft(channel_est_time(1,:));
            
            % FFT Preamble
            preamble_out_freq1 = fft(preamble_out1);
            preamble_out_freq2 = fft(preamble_out2);
            
            channel_est1 = preamble_out_freq1./preamble;
            channel_est2 = preamble_out_freq2./preamble;
            channel_est = (channel_est1 + channel_est2)/2;
           
            
            if estimation == true
                if c == 3 
                    channels{3} = channel_est_timefreq;
                end
                if c == 4
                    channels{4} = channel_est;
                end
            end       
                
            % S/P
            y_rec = reshape(y_noisy_nopre, [N + cp_length, length(y_noisy_nopre)/(N + cp_length)]);
            
            % CP Removal
            signal_no_cp = y_rec(cp_length + 1:end, :);

            % FFT
            freq_signal = fft(signal_no_cp);

            % Equalization
            equalized_signal = zeros(size(freq_signal));
            if c == 1 || c == 2
                for j = 1:size(freq_signal, 2)
                    equalized_signal(:,j) = freq_signal(:,j) ./ fft(channels{c}, size(freq_signal, 1)).'; % True channel
                end
            else
                for j = 1:size(freq_signal, 2)
                    equalized_signal(:,j) = freq_signal(:,j) ./ channels{c}.'; % Channel esimation
                end
            end
            equalized_signal = zeros(size(freq_signal));
            for j = 1:size(freq_signal, 2)
                equalized_signal(:,j) = freq_signal(:,j) ./ channel_est.'; % Channel esimation
            end
            % P/S
            serial_signal = reshape(equalized_signal, [1, numel(equalized_signal)]);

            % Demapping
            bit_rx = demapping1(serial_signal.',Nbps,modulation);

            % BER 
            diffs = abs(bit_tx' - bit_rx);
            errors = sum(diffs);
            BER = errors/length(bit_tx);
            sumBER = sumBER + BER;
            
            % NMSE
            if c == 3 
                NMSE = norm(channel_est_timefreq - fft(channels{1},size(channel_est_timefreq,2)),2)/norm(fft(channels{1},size(channel_est_timefreq,2)),2);
            elseif c == 4
                NMSE = norm(channel_est - fft(channels{2},size(channel_est,2)),2)/norm(fft(channels{2},size(channel_est,2)),2);
            else
                NMSE = 0;
            end
            sumNMSE = sumNMSE + NMSE;
            
            % RMSE
            RMSE = (Toa - toa_est)^2;
            sumRMSE = sumRMSE + RMSE;
            
        end
        ber_vec(i,c) = sumBER/essais;
        nmse_vec(i,c) = sumNMSE/essais;
        rmse_vec(i,c) = sqrt(sumRMSE/essais);
        
    end
end

%% Plot BER
figure
% semilogy(Eb_N0_dB,ber_vec(:,1),'r');
% hold on;
% semilogy(Eb_N0_dB,ber_vec(:,2),'b');
% semilogy(Eb_N0_dB,ber_vec(:,3),'r--');
% semilogy(Eb_N0_dB,ber_vec(:,4),'b--');
% grid on
% 
% xlabel('Eb/N0 (dB)')
% ylabel('BER')
% grid on;
% legend('NLOS no TOA - with channel estimation','LOS no TOA - with channel estimation','NLOS with TOA and channel estimation','LOS with TOA and channel estimation')

%% Plot RMSE
% figure;
% plot(Eb_N0_dB,rmse_vec(:,1),'b');
% grid on;
% xlabel('Eb/N0 (dB)');
% ylabel('RMSE of the TOA estimate');
% save('cross','rmse_vec');