%% LABO 6

clear all;
close all;
clc;

% Load impulse responses
data = load('Stat_channel.mat');
h_wide_NLOS_stat = data.h_NLOS_stat;
h_wide_LOS_stat = data.h_LOS_stat;

%% PARAMETERS
avgtot = 0;
Nbps = 4;
nb_bits = 4*3840;
bit_tx = randi([0,1],[1 nb_bits]);
Fs = 400e3;
Ts = 1/Fs;
modulation = 'qam';
N = 64;
cp_length = 16;
Time_acquisition = false;
Frequency_acquisition = true;
phase_tracking = false;
%CFO = randi([1 100])/200;
CFO = 1000;
%phase_error = randi([0 1000])*pi/1000;
phase_error = 0;
first_time = 0;
CFOvec = [0 1e3 2e3 3e3];

%% TRANSMITTER
symb = mapping(bit_tx',Nbps,modulation);
%scatterplot(symb)

% S/P
stream = reshape(symb, [N length(symb)/(N)]);

if Frequency_acquisition == true
    stream = reshape(symb, [N-4 length(symb)/(N-4)]);
    stream1 = zeros(64,64);
    for i = 1:64
        str = stream(:,i).';
        res = [str(1:10),1,str(11:23),1,str(24:36),1,str(37:49),1,str(50:end)];
        stream1(:,i) = res.'; 
    end
    stream = stream1;
end

%plot(abs(res))

% IFFT
stream_time = ifft(stream);

% CP
if Frequency_acquisition == 1
    symb_ifft_cp = zeros(N + cp_length, size(stream_time,1));
else 
    symb_ifft_cp = zeros(N + cp_length, length(symb)/N);
end
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
if Time_acquisition == true
    Toa = randi(round(length(frame)*0.1));
else
    Toa = 0;
end


%% CHANNEL
vec = [false true];
channel = h_wide_NLOS_stat;
Eb_N0_dB = -5:25;
ber_vec = zeros(length(vec),length(Eb_N0_dB));
nmse_vec = zeros(1,length(Eb_N0_dB));
essais = 20;
% CFO_vec = [0 1e3 2e3 3e3];

    
% Convolution with Channel
y = conv(channel,frame);
y = y(1:end-5);
if Time_acquisition == true
    y = [zeros(1,Toa) y zeros(1,2000)];
end

% for cc = 1:length(vec)
for cc = 1:length(CFOvec)
    cc
    CFO = CFOvec(cc)
    %phase_tracking = vec(cc);
    for i = 1:length(Eb_N0_dB)
        Eb_N0_dB(i);
        sumBER = 0;
        sumNMSE = 0;
        for tries = 1:essais
            %% AWG Noise
            Eb = var(y)/Nbps;
            EbN0 =  10^(Eb_N0_dB(i)/10); 
            N0 = Eb/EbN0;
            noise = sqrt(N0/2)*(randn(1,length(y))+1i*randn(1,length(y)));
            y_noisy = y + noise;

            %% Add CFO
            if Frequency_acquisition == true
                axis = 0:Ts:Ts*(length(y_noisy)-1);
                y_noisy = y_noisy.*exp(1j*(2*pi*CFO*axis+phase_error));
            end

            %% RECEIVER
            lengthCorr = round(length(frame)*0.1);
            correlation_vector = zeros(1,lengthCorr);
            somme = zeros(1,lengthCorr);

            for n = 1:lengthCorr
                for l = 0:N-1
                    correlation_vector(n) =  correlation_vector(n) + y_noisy(n+l+64).* conj(y_noisy(n+l));
                end
                for l = 0: 2*N -1
                    somme(n) = somme (n) + abs(y_noisy(n +l)).^2;
                end
            end

            CFO_est = (phase(correlation_vector));

            correlation = abs(correlation_vector)./(somme);
    %         if i == 10
    %             figure;
    %             plot(correlation);
    %             title('Correlation');
    %         end
            M = zeros(1,length(correlation));
            AvgSize = 32;

            for k = 1:length(correlation)-AvgSize
                M(k) = sum(abs(correlation(k:k+AvgSize)))/(AvgSize+1);
            end

    %         if i == 10
    %             figure;
    %             plot(M);
    %             title('Correlation with Moving Average Window');
    %         end

            [ ~, delayEst] = max (M);
            CFO_est = CFO_est(delayEst+cp_length);

            if CFO_est<0
                CFO_est = CFO_est + 2*pi;
            end
            if CFO_est>pi
                CFO_est = CFO_est - 2*pi;
            end

            CFO_est = CFO_est*1000;

    %         toa_est = delayEst-15;
    % 
    %         if(toa_est<=0)
    %            toa_est = 1; 
    %         end
    %         y_noisy = y_noisy(toa_est : toa_est + length(frame)-1);     

            % Frequency acquisition
            if Frequency_acquisition == true
                axis = 0:Ts:Ts*(length(y_noisy)-1);
                y_noisy = y_noisy.*exp(-1j*(2*pi*abs(CFO_est)*axis));
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


            % S/P
            y_rec = reshape(y_noisy_nopre, [N + cp_length, length(y_noisy_nopre)/(N + cp_length)]);

            % CP Removal
            signal_no_cp = y_rec(cp_length + 1:end, :);

            % FFT
            freq_signal = fft(signal_no_cp);

            % Equalization
            equalized_signal = zeros(size(freq_signal));
            for j = 1:size(freq_signal, 2)
                equalized_signal(:,j) = freq_signal(:,j) ./ channel_est.'; % Channel esimation
            end
            avgtot = 0;
            if phase_tracking == true
                for m = 1:size(equalized_signal,2)
                    %Tracking
                    a = phase(equalized_signal(53,m))- phase(equalized_signal(39,m));
                    b = phase(equalized_signal(39,m))- phase(equalized_signal(25,m));
                    x = phase(equalized_signal(25,m))- phase(equalized_signal(11,m));
    %                     d = phase(equalized_signal(11,m))- phase(equalized_signal(53,m));

                    phase_err_awg = ((a+b+x)/3);
                    avgtot = avgtot + phase_err_awg;
    %                         equalized_signal(:,m) = equalized_signal(:,m).*exp(-1j*phase_err_awg);
                end
                avgtot = avgtot/size(equalized_signal,2);
                equalized_signal = equalized_signal.*exp(-1j*2*pi*avgtot);
            end

           if Frequency_acquisition == true
                equalized_signal1 = zeros(N-4,N);
                for p = 1:N
                    str = equalized_signal(:,p).';
                    res = [str(1:10),str(12:24),str(26:38),str(40:52),str(54:end)];
                    equalized_signal1(:,p) = res.'; 
                end
                equalized_signal = equalized_signal1;
           end 

            % P/S
            serial_signal = reshape(equalized_signal, [1, numel(equalized_signal)]);

    %         if Frequency_acquisition == true
    %             if first_time == 0
    %                 scatterplot(serial_signal);
    %                 first_time = first_time+1;
    %             end
    %         end

            % Demapping
            bit_rx = demapping1(serial_signal.',Nbps,modulation);

            % BER 
            diffs = abs(bit_tx' - bit_rx);
            errors = sum(diffs);
            BER = errors/length(bit_tx);
            sumBER = sumBER + BER;

        end

        ber_vec(cc,i) = sumBER/essais;
    end

end

%% Plot BER CFO
figure
semilogy(Eb_N0_dB,ber_vec(1,:),'.-');hold on;
semilogy(Eb_N0_dB,ber_vec(2,:),'.-');hold on;
semilogy(Eb_N0_dB,ber_vec(3,:),'.-');hold on;
semilogy(Eb_N0_dB,ber_vec(4,:),'.-');hold on;
grid on

xlabel('Eb/N0 (dB)')
ylabel('BER')
grid on;
title('NLOS channel, 16 QAM')
legend('No CFO','CFO = 1000 Hz','CFO = 2000 Hz','CFO = 3000 Hz');
%legend('NLOS with TOA,CFO and channel estimation','LOS with TOA,CFO and channel estimation','NLOS with only channel estimation','LOS with only channel estimation')
%save('ber_cfo300','ber_vec');

%% Plot BER phase tracking
% figure
% semilogy(Eb_N0_dB,ber_vec(1,:),'.-');hold on;
% semilogy(Eb_N0_dB,ber_vec(2,:),'.-');hold on;
% xlabel('Eb/N0 (dB)')
% ylabel('BER')
% grid on;
% title('Phase tracking of CFO, NLOS channel, 16 QAM')
% legend('Acquisition only','Acquisition + Phase tracking');
% grid on
%%
% a = load('ber_cfo10');
% b = load('ber_cfo20');
% c = load('ber_cfo30');
% a = a.ber_vec;
% b = b.ber_vec;
% c = c.ber_vec;
% figure
% semilogy(Eb_N0_dB,a,'.-');hold on;
% semilogy(Eb_N0_dB,b,'.-');hold on;
% semilogy(Eb_N0_dB,c,'.-');hold on;
% grid on
% xlabel('Eb/N0 (dB)')
% ylabel('BER')
% grid on;
% title('NLOS channel, 16 QAM')
% legend('CFO = 10','CFO = 20','CFO = 30')

