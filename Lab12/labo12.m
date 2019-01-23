%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project : Conception of a complete OFDM             % 
% Authors : Alexis Rouvroy, Léo Moulin & Rudy Schoone %
% Practicals sessions 1 & 2                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data files

close all
clear all 
clc

LOS_data = load('LOS.mat');
NLOS_data = load('NLOS.mat');

x = LOS_data.xx;
y = LOS_data.yy;
z = LOS_data.zz;

%% Initialize data

Len_data_wideband = length(NLOS_data.Data{x}{y}{z});   % Len_data_wideband   = 501
Len_data_narrowband = ceil(Len_data_wideband/10);      % Len_data_narrowband = 51
Ts = 1/200e6;       % Sampling period : fs = 200 MHz
PDP = 50;
PDP_down = 20;
PDP_NLOS = 0;
PDP_LOS = 0;
PDP_NLOS_down = 0;
PDP_LOS_down = 0;
i=1;

% Create rectangular window

w = rectwin(Len_data_narrowband)';  % w = ones(51,1)
w = padarray(w,[0 225]);            % add 225 zeros before and after the windiw)
% figure;
% plot(w);
% title('Rectangular window');

%% Compute data
% Initialyze the matrix
H_NLOS      = zeros(x,y,z,Len_data_wideband);
H_LOS       = zeros(x,y,z,Len_data_wideband);
H_NLOS_down = zeros(x,y,z,Len_data_wideband);
H_LOS_down  = zeros(x,y,z,Len_data_wideband);

h_save_LOS  = zeros(51,1000);
h_save_NLOS = zeros(51,1000);

h_narrow_NLOS = zeros(1000,1);
h_narrow_LOS  = zeros(1000,1);

h_wide_NLOS = zeros(1000,10);
h_wide_LOS  = zeros(1000,10);

for xx=1:x
    for yy=1:y
        for zz =1:z
            % frequency response
            H_NLOS(xx,yy,zz,:)      = NLOS_data.Data{xx}{yy}{zz};
            H_LOS(xx,yy,zz,:)       = LOS_data.Data{xx}{yy}{zz};
            H_NLOS_down(xx,yy,zz,:) = (NLOS_data.Data{xx}{yy}{zz}).*w; % Apply the rectangular window on the data
            H_LOS_down(xx,yy,zz,:)  = (LOS_data.Data{xx}{yy}{zz}).*w;

            % Impulse response
            h_NLOS      = sqrt(length(H_NLOS))*ifft(ifftshift(H_NLOS(xx,yy,zz,:)));
            h_LOS       = sqrt(length(H_LOS))*ifft(ifftshift(H_LOS(xx,yy,zz,:)));    
            h_NLOS_down = sqrt(length(H_NLOS_down))*ifft(ifftshift(H_NLOS_down(xx,yy,zz,:)));
            h_LOS_down  = sqrt(length(H_LOS_down))*ifft(ifftshift(H_LOS_down(xx,yy,zz,:)));              
            
            % Downsample to get the narrowband signal
            h_NLOS_down = downsample(h_NLOS_down, 10);
            h_LOS_down  = downsample(h_LOS_down, 10);
            
            % Compute indice of maximum value
            [~,I_NLOS]      = max(h_NLOS(1:PDP));
            [~,I_LOS]       = max(h_LOS(1:PDP)); 
            [~,I_NLOS_down] = max(h_NLOS_down(1:PDP_down));
            [~,I_LOS_down]  = max(h_LOS_down(1:PDP_down));  

            % 4D to 1D matrix (vector)
            h_NLOS           = h_NLOS(:);
            h_LOS            = h_LOS(:);
            h_LOS_down       = h_LOS_down(:);
            h_NLOS_down      = h_NLOS_down(:);
            h_save_LOS(:,i)  = h_LOS_down;
            h_save_NLOS(:,i) = h_NLOS_down;

            % Compute PDP
            PDP_NLOS      = abs(h_NLOS(I_NLOS:I_NLOS+PDP)).^2 + PDP_NLOS;          % Sum the PDP of each measurement     
            PDP_LOS       = abs(h_LOS(I_LOS:I_LOS+PDP)).^2 + PDP_LOS;            
            PDP_NLOS_down = abs(h_NLOS_down(I_NLOS_down:I_NLOS_down+PDP_down)).^2 + PDP_NLOS_down;            
            PDP_LOS_down  = abs(h_LOS_down(I_LOS_down:I_LOS_down+PDP_down)).^2 + PDP_LOS_down;         
            
            % Narrowband
            h_narrow_NLOS(i,:) = abs(sum(h_NLOS_down(1: 1+6)));
            h_narrow_LOS(i,:)  = abs(sum(h_LOS_down(1: 1+6)));

            % Wideband
            h_wide_NLOS(i,:) = abs(h_NLOS_down(1:10));
            h_wide_LOS(i,:)  = abs(h_LOS_down(1:10));
            
            i = i+1;
        end
    end
end

PDP_LOS       = PDP_LOS/1000;    % "Normalize" the PDP
PDP_NLOS      = PDP_NLOS/1000;
PDP_NLOS_down = PDP_NLOS_down/1000;
PDP_LOS_down  = PDP_LOS_down/1000;

% Compute Delay spread
mean_delay_NLOS = 0;
mean_delay_LOS  = 0;
tau_NLOS = 0;
tau_LOS  = 0;

for i = 1:PDP + 1
    % NLOS
    mean_delay_NLOS = mean_delay_NLOS + i*Ts*PDP_NLOS(i);
    tau_NLOS = tau_NLOS + (i*Ts)^2 * PDP_NLOS(i);
    
    % LOS
    mean_delay_LOS = mean_delay_LOS + i*Ts*PDP_LOS(i);
    tau_LOS = tau_LOS + (i*Ts)^2 * PDP_LOS(i);
end

mean_delay_NLOS_down = 0;
mean_delay_LOS_down  = 0;
tau_NLOS_down = 0;
tau_LOS_down  = 0;

for i = 1:PDP_down + 1
    % NLOS
    mean_delay_NLOS_down = mean_delay_NLOS_down + i*Ts*PDP_NLOS_down(i);
    tau_NLOS_down = tau_NLOS_down + (i*Ts)^2 * PDP_NLOS_down(i);
    
    % LOS
    mean_delay_LOS_down = mean_delay_LOS_down + i*Ts*PDP_LOS_down(i);
    tau_LOS_down = tau_LOS_down + (i*Ts)^2 * PDP_LOS_down(i);
end

% LOS
total_power_LOS = sum(PDP_LOS);
mean_delay_LOS  = mean_delay_LOS/total_power_LOS;
sigma_tau_LOS   = sqrt(tau_LOS/total_power_LOS - mean_delay_LOS^2); % Delay spread

% LOS reduced
total_power_LOS_down = sum(PDP_LOS_down);
mean_delay_LOS_down  = mean_delay_LOS_down/total_power_LOS_down;
sigma_tau_LOS_down   = sqrt(tau_LOS_down/total_power_LOS_down - mean_delay_LOS_down^2); % Delay spread

% NLOS
total_power_NLOS = sum(PDP_NLOS);
mean_delay_NLOS  = mean_delay_NLOS/total_power_NLOS;
sigma_tau_NLOS   = sqrt(tau_NLOS/total_power_NLOS - mean_delay_NLOS^2); % Delay spread

% NLOS reduced
total_power_NLOS_down = sum(PDP_NLOS_down);
mean_delay_NLOS_down  = mean_delay_NLOS_down/total_power_NLOS_down;
sigma_tau_NLOS_down   = sqrt(tau_NLOS_down/total_power_NLOS_down - mean_delay_NLOS_down^2); % Delay spread

% PDP_LOS modeled by an exponential
tau_LOS = (0:PDP) * Ts;
P_tau_exp_LOS = PDP_LOS(1) * exp(-tau_LOS/sigma_tau_LOS);

% PDP_LOS_down modeled by an exponential
tau_LOS_down = (0:PDP_down) * Ts;
P_tau_exp_LOS_down = PDP_LOS_down(1) * exp(-tau_LOS_down/sigma_tau_LOS_down);

% PDP_NLOS modeled by an exponential
tau_NLOS = (0:PDP) * Ts;
P_tau_exp_NLOS = PDP_NLOS(1) * exp(-tau_NLOS/sigma_tau_NLOS);

% PDP_NLOS_down modeled by an exponential
tau_NLOS_down = (0:PDP_down) * Ts;
P_tau_exp_NLOS_down = PDP_NLOS_down(1) * exp(-tau_NLOS_down/sigma_tau_NLOS_down);

% Compute coherence bandwith
Delta_fc_NLOS = 1/(2*pi*sigma_tau_NLOS);
Delta_fc_LOS  = 1/(2*pi*sigma_tau_LOS);
Delta_fc_NLOS_down = 1/(2*pi*sigma_tau_NLOS_down);
Delta_fc_LOS_down  = 1/(2*pi*sigma_tau_LOS_down);

%% Print figures

figure
subplot(2,2,1)
plot(10*log10(PDP_NLOS))
hold on
plot(10*log10(P_tau_exp_NLOS),'r')
title('PDP-NLOS')
legend('PDP','PDP exp')
xlabel('Time (s)'); 

subplot(2,2,2)
plot(10*log10(PDP_LOS))
hold on
plot(10*log10(P_tau_exp_LOS),'r')
title('PDP-LOS')
legend('PDP','PDP exp')
xlabel('Time (s)'); 

subplot(2,2,3)
plot(10*log10(PDP_LOS_down))
hold on
plot(10*log10(P_tau_exp_LOS_down),'r')
title('PDP-LOS-down')
legend('PDP','PDP exp')
xlabel('Time (s)'); 

subplot(2,2,4)
plot(10*log10(PDP_NLOS_down))
hold on
plot(10*log10(P_tau_exp_NLOS_down),'r')
title('PDP-NLOS-down')
legend('PDP','PDP exp')
xlabel('Time (s)'); 

% Display results
disp(['Delay spread_LOS = ', num2str(sigma_tau_LOS)])
disp(['Coherence bandwidth_LOS = ', num2str(Delta_fc_LOS)])
disp(['Delay spread_NLOS = ', num2str(sigma_tau_NLOS)])
disp(['Coherence bandwidth_NLOS = ', num2str(Delta_fc_NLOS)])
disp(['Delay spread_LOS_down = ', num2str(sigma_tau_LOS_down)])
disp(['Coherence bandwidth_LOS_down = ', num2str(Delta_fc_LOS_down)])
disp(['Delay spread_NLOS_down = ', num2str(sigma_tau_NLOS_down)])
disp(['Coherence bandwidth_NLOS_down = ', num2str(Delta_fc_NLOS_down)])

%% LAB 2

%% Narrowband model

% Narrowband NLOS

% dfittool(h_narrow_NLOS)
% s_NLOS_narrow     = 1.52286e-06;
s_NLOS_narrow     = 1.59427e-06;
% sigma_NLOS_narrow = 0.00193872;
sigma_NLOS_narrow = 0.0019273;

K_NLOS_narrow = 10*log10((s_NLOS_narrow)^2/(2*sigma_NLOS_narrow^2)); % No LOS => K = 0

% Narrowband Statistical model 
distr_NLOS_wide    = makedist('Rician','s',s_NLOS_narrow,'sigma',sigma_NLOS_narrow); % create a Rician (Rayleigh) distribution 
A_narrow_NLOS      = random(distr_NLOS_wide);  
phi_narrow_NLOS    = rand(1)*2*pi;
h_narrow_NLOS_stat = A_narrow_NLOS*exp(1j*phi_narrow_NLOS);

disp(' ')
disp('Narrowband NLOS:')
disp(['K_NLOS_narrow = ', num2str(K_NLOS_narrow)])
disp(['Statistical Model Narrowband NLOS = ', num2str(h_narrow_NLOS_stat)])

% Narrowband LOS

% dfittool(h_narrow_LOS)
s_LOS_narrow     = 0.0128789;
sigma_LOS_narrow = 0.00544533;

K_LOS_narrow = 10*log10((s_LOS_narrow)^2/(2*sigma_LOS_narrow^2));

% Narrowband Statistical model 
distr_LOS         = makedist('Rician','s',s_LOS_narrow,'sigma',sigma_LOS_narrow);
A_narrow_LOS      = random(distr_LOS);
phi_narrow_LOS    = rand(1)*2*pi;
h_narrow_LOS_stat = A_narrow_LOS*exp(1j*phi_narrow_LOS);

disp(' ')
disp('Narrowband LOS:')
disp(['K_LOS_narrow = ', num2str(K_LOS_narrow)])
disp(['Statistical Model Narrowband LOS = ', num2str(h_narrow_LOS_stat)])

%% Wideband model
% Matrix initialization
K_NLOS_wide      = zeros(1,length(h_wide_NLOS(1,:))); 
h_wide_NLOS_stat = zeros(1,length(h_wide_NLOS(1,:)));
K_LOS_wide       = zeros(1,length(h_wide_NLOS(1,:))); 
h_wide_LOS_stat  = zeros(1,length(h_wide_NLOS(1,:)));

for i = 1:length(h_wide_NLOS(1,:))
    % NLOS
    Wideband_NLOS   = fitdist(h_wide_NLOS(:,i),'Rician');
    s_NLOS_wide     = Wideband_NLOS.s;
    sigma_NLOS_wide = Wideband_NLOS.sigma;
    K_NLOS_wide(i)  = 10*log10((s_NLOS_wide)^2/(2*sigma_NLOS_wide^2));
    % Statistichal distribution
    distr_NLOS_wide     = makedist('Rician','s',s_NLOS_wide,'sigma',sigma_NLOS_wide);
    A_wide_NLOS         = random(distr_NLOS_wide);
    phi_wide_NLOS       = rand(1)*2*pi;
    h_wide_NLOS_stat(i) = A_wide_NLOS*exp(1j*phi_wide_NLOS);
    
    % LOS
    Wideband_LOS   = fitdist(h_wide_LOS(:,i),'Rician');
    s_LOS_wide     = Wideband_LOS.s;
    sigma_LOS_wide = Wideband_LOS.sigma;
    K_LOS_wide(i)  = 10*log10((s_LOS_wide)^2/(2*sigma_LOS_wide^2));
    % Statistichal distribution
    distr_LOS_wide     = makedist('Rician','s',s_LOS_wide,'sigma',sigma_LOS_wide);
    A_wide_LOS         = random(distr_LOS_wide);
    phi_wide_LOS       = rand(1)*2*pi;
    h_wide_LOS_stat(i) = A_wide_LOS*exp(1j*phi_wide_LOS);
end

%% Print figures

% NLOS
figure
plot(K_NLOS_wide)
title('Evolution of K as a function of the delay-NLOS')
ylabel('K (dB)');
xlabel('Delay (tap)'); 

% Wideband Statistical model --> sum of the narrowband model
disp(['Statistical Model Wideband = ', num2str(h_wide_NLOS_stat)])
figure
plot(abs(h_wide_NLOS_stat))
title('Statistical model Wideband in function of taps-NLOS');

% LOS

figure
plot(K_LOS_wide)
title('Evolution of K as a function of the delay-LOS')
ylabel('K (dB)');
xlabel('Delay (tap)'); 

% Wideband Statistical model --> sum of the narrowband model
disp(['Statistical Model Wideband = ', num2str(h_wide_LOS_stat)])
figure
stem(abs(h_wide_LOS_stat))
hold on
stem(abs(h_wide_LOS(2,:)),'r')
legend('Statistical model','Measurements');
title('Statistical model Wideband in function of taps-LOS');

h_save_LOS = h_save_LOS(1:6,:).';
h_save_NLOS = h_save_NLOS(2:7,:).';
save('test.mat','h_save_LOS','h_save_NLOS'); % Saving paramters

h_NLOS_stat = h_wide_NLOS_stat(2:7);
h_LOS_stat = h_wide_LOS_stat(1:6);
save('Stat_channel.mat','h_LOS_stat','h_NLOS_stat');

figure
plot(abs(h_NLOS_stat))
hold on
plot(abs(h_LOS_stat),'r')







