%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      LAB 7     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    Beamforming method     %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Moulin Léo - Rouvroy Alexis - Schoone Rudy  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all

%% Initialize parameters

% Axis position of the different antennas
X=0:9;  
Z=0:9;
Y=9:-1:0;
% Bandwidth frequency
freq=2.35e9;    
lambda=3e8/freq;
% Angle for the receiving waves
step=0.1;
phi=-pi:step:pi;
theta=0:step:pi;
% Load old channels
h_old = load('Channels_old.mat');

h_LOS_narrow = h_old.h_LOS_narrow;
h_NLOS_narrow = h_old.h_NLOS_narrow;


%% 3D Beamforming function

[A_LOS,A_NLOS] = Beamforming(X,Y,Z,theta,phi,lambda,h_LOS_narrow,h_NLOS_narrow);

%% Plot LOS case
angularSpectrum=20*log10(reshape(abs(A_LOS(1,:,:)),[length(theta) length(phi)]));
Sup=max(max(angularSpectrum));
Inf=min(min(angularSpectrum));

% Plot the angular spectrum on the image of the Lab
plotImage(phi,theta,angularSpectrum,[Sup Inf])
colorbar('YTickLabel',{'-110','-90','-70','-50'});
colormap jet;hold on;
title('Angular spetrum in the LOS case for the 20MHz channel');
xlabel('Phi');ylabel('Theta');

% Find the local Maxima of the angular spectrum
index_maximum = findLocalMaxima(angularSpectrum,-50);
plot(phi(index_maximum(2,:)),theta(index_maximum(1,:)),'*');
figure;

% Plot the angular spectrum in 3D
surf(reshape(abs(A_LOS(1,:,:)),[length(theta) length(phi)]));
xlabel('Phi');ylabel('Theta');zlabel('|a(n)|');

%% Plot NLOS case
angularSpectrum=20*log10(reshape(abs(A_NLOS(2,:,:)),[length(theta) length(phi)]));
Sup=max(max(angularSpectrum));
Inf=min(min(angularSpectrum));

% Plot the angular spectrum on the image of the Lab
plotImage(phi,theta,angularSpectrum,[Sup Inf])
colorbar('YTickLabel',{'-130','-110','-90','-70'});
hold on;
title('Angular spetrum in the NLOS case for the 20MHz channel');
xlabel('Phi');ylabel('Theta');

% Find the local Maxima of the angular spectrum
index_maximum = findLocalMaxima(angularSpectrum,-70);
plot(phi(index_maximum(2,:)),theta(index_maximum(1,:)),'*');
figure;

% Plot the angular spectrum in 3D
surf(reshape(abs(A_NLOS(5,:,:)),[length(theta) length(phi)]));
xlabel('Phi');ylabel('Theta');zlabel('|a(n)|');

%% Build a new channel 

%  [H_new_LOS, H_new_NLOS] = Compute_Channel(A_LOS,A_NLOS,theta,phi,X,Y,Z,lambda);
%  save('Channels_new.mat','H_new_LOS','H_new_NLOS'); % Saving paramters
%% Load new channel (Computation of new channel takes time)

H_new = load('Channels_new.mat');
H_new_LOS = H_new.H_new_LOS;
H_new_NLOS = H_new.H_new_NLOS;
H_new_LOS_2D = reshape(abs(H_new_LOS(1,1,1,1:6)),[1 6]);
H_new_NLOS_2D = reshape(abs(H_new_NLOS(1,1,1,2:7)),[1 6]);
figure
plot(H_new_LOS_2D)
hold on
plot(H_new_NLOS_2D,'r')
xlabel('tap');ylabel('|h(n)|');
title('Computation of the new channel model based on the beamforming')
legend('LOS case','NLOS case')

%%
figure
plot(abs(h_LOS_narrow(1,1:6)))
hold on 
plot(abs(h_NLOS_narrow(1,2:7)),'r')
legend('LOS case', 'NLOS case')
title('Old channel model')

%% Spatial correlation
deltaX = 0:0.0002:0.2;
deltaY = 0:0.0002:0.2;
deltaZ = [0:9]*0.02;
beta = 2*pi/lambda;
deltaZ = 0:0.0002:0.2;
R_deltaZ_LOS = zeros(length(H_new_LOS(1,1,1,:)),length(deltaZ));
R_deltaZ_NLOS = zeros(length(H_new_NLOS(1,1,1,:)),length(deltaZ));
n=1;
tap=1;
threshold = 1e-4;
newAngles_LOS = findLocalMaxima(abs(squeeze(A_LOS(tap,:,:))), threshold);
newAngles_NLOS = findLocalMaxima(abs(squeeze(A_NLOS(tap,:,:))), threshold);

a_LOS = squeeze(A_LOS(tap,:,:));
a_NLOS = squeeze(A_NLOS(tap,:,:));
a_theta_LOS = mean(a_LOS(:,:),2);
a_theta_NLOS = mean(a_NLOS(:,:),2);
RZ_LOS = zeros(size(deltaZ));
RZ_NLOS = zeros(size(deltaZ));
for z = 1:length(deltaZ)
    RZ_LOS(z) = mean(abs(a_theta_LOS(:)).^2 .* exp(1i*beta.*cos(theta).*deltaZ(z))' );
    RZ_NLOS(z) = mean(abs(a_theta_NLOS(:)).^2 .* exp(1i*beta.*cos(theta).*deltaZ(z))' );
end
RX_LOS = zeros(size(deltaX));
RX_NLOS = zeros(size(deltaX));
for x = 1:length(deltaX)
    RX_LOS(x) = mean(abs(a_theta_LOS(newAngles_LOS(1,:))).*exp(1i*beta.*sin(newAngles_LOS(1,:)).*cos(newAngles_LOS(2,:)).*deltaX(x)).');
    RX_NLOS(x) = mean(abs(a_theta_NLOS(newAngles_NLOS(1,:))).*exp(1i*beta.*sin(newAngles_NLOS(1,:)).*cos(newAngles_NLOS(2,:)).*deltaX(x)).');
end
RY_LOS = zeros(size(deltaY));
RY_NLOS = zeros(size(deltaY));
for y = 1:length(deltaY)
    RY_LOS(y) = mean(abs(a_theta_LOS(newAngles_LOS(1,:))) .* exp(1i*beta.*sin(newAngles_LOS(1,:)).*sin(newAngles_LOS(2,:)).*deltaY(y)).');
    RY_NLOS(y) = mean(abs(a_theta_NLOS(newAngles_NLOS(1,:))) .* exp(1i*beta.*sin(newAngles_NLOS(1,:)).*sin(newAngles_NLOS(2,:)).*deltaY(y)).');
end

RZ_LOS=abs(RZ_LOS)/abs(max(RZ_LOS));
RX_LOS=abs(RX_LOS)/abs(max(RX_LOS));
RY_LOS=abs(RY_LOS)/abs(max(RY_LOS));
RZ_NLOS=abs(RZ_NLOS)/abs(max(RZ_NLOS));
RX_NLOS=abs(RX_NLOS)/abs(max(RX_NLOS));
RY_NLOS=abs(RY_NLOS)/abs(max(RY_NLOS));

rtot=RX_LOS+RY_LOS+RZ_LOS;
rtot=rtot/max(rtot);
a = (abs(besselj(0,2*pi/lambda*deltaZ)));
figure
plot(deltaX/lambda,RX_LOS,'b');
hold on 
plot(deltaY/lambda,RY_LOS,'g');
hold on 
plot(deltaZ/lambda,RZ_LOS,'r');
grid on
plot(deltaY/lambda,a,'m');
xlabel('distance[m]');ylabel('Correlation');
title('Spatial correlation of the new channel model for each direction in wideband - LOS case')
legend('X LOS case','Y LOS case','Z LOS case','clarke model')

figure
plot(deltaX/lambda,RX_NLOS,'b');
hold on 
plot(deltaY/lambda,RY_NLOS,'g');
hold on 
plot(deltaZ/lambda,RZ_NLOS,'r');
grid on
hold on
plot(deltaY/lambda,a,'m');
xlabel('distance[m]');ylabel('Correlation');
title('Spatial correlation of the new channel model for each direction in wideband - NLOS case')
legend('X NLOS case','Y NLOS case','Z NLOS case','clarke model')



%% Narrowband
A_save_LOS = 0;
A_save_NLOS = 0;
for i = 1:6
    A_save_LOS = A_save_LOS + A_LOS(i,:,:);
    A_save_NLOS = A_save_NLOS + A_NLOS(i,:,:);
end

a_LOS = squeeze(A_save_LOS(tap,:,:));
a_NLOS = squeeze(A_save_NLOS(tap,:,:));
a_theta_LOS = mean(a_LOS(:,:),2);
a_theta_NLOS = mean(a_NLOS(:,:),2);
RZ_LOS = zeros(size(deltaZ));
RZ_NLOS = zeros(size(deltaZ));
for z = 1:length(deltaZ)
    RZ_LOS(z) = mean(abs(a_theta_LOS(:)).^2 .* exp(1i*beta.*cos(theta).*deltaZ(z))' );
    RZ_NLOS(z) = mean(abs(a_theta_NLOS(:)).^2 .* exp(1i*beta.*cos(theta).*deltaZ(z))' );
end
RX_LOS = zeros(size(deltaX));
RX_NLOS = zeros(size(deltaX));
for x = 1:length(deltaX)
    RX_LOS(x) = mean(abs(a_theta_LOS(newAngles_LOS(1,:))).*exp(1i*beta.*sin(newAngles_LOS(1,:)).*cos(newAngles_LOS(2,:)).*deltaX(x)).');
    RX_NLOS(x) = mean(abs(a_theta_NLOS(newAngles_NLOS(1,:))).*exp(1i*beta.*sin(newAngles_NLOS(1,:)).*cos(newAngles_NLOS(2,:)).*deltaX(x)).');
end
RY_LOS = zeros(size(deltaY));
RY_NLOS = zeros(size(deltaY));
for y = 1:length(deltaY)
    RY_LOS(y) = mean(abs(a_theta_LOS(newAngles_LOS(1,:))) .* exp(1i*beta.*sin(newAngles_LOS(1,:)).*sin(newAngles_LOS(2,:)).*deltaY(y)).');
    RY_NLOS(y) = mean(abs(a_theta_NLOS(newAngles_NLOS(1,:))) .* exp(1i*beta.*sin(newAngles_NLOS(1,:)).*sin(newAngles_NLOS(2,:)).*deltaY(y)).');
end
RZ_LOS=abs(RZ_LOS)/abs(max(RZ_LOS));
RX_LOS=abs(RX_LOS)/abs(max(RX_LOS));
RY_LOS=abs(RY_LOS)/abs(max(RY_LOS));
RZ_NLOS=abs(RZ_NLOS)/abs(max(RZ_NLOS));
RX_NLOS=abs(RX_NLOS)/abs(max(RX_NLOS));
RY_NLOS=abs(RY_NLOS)/abs(max(RY_NLOS));

a = (abs(besselj(0,2*pi/lambda*deltaY)));
figure
plot(deltaX/lambda,RX_LOS,'b');
hold on 
plot(deltaY/lambda,RY_LOS,'g');
hold on 
plot(deltaZ/lambda,RZ_LOS,'r');
grid on
plot(deltaY/lambda,a,'m');
xlabel('distance[m]');ylabel('Correlation');
title('Spatial correlation of the new channel model for each direction in narrowband - LOS case')
legend('X LOS case','Y LOS case','Z LOS case','clarke model')


figure
plot(deltaX/lambda,RX_NLOS,'b');
hold on 
plot(deltaY/lambda,RY_NLOS,'g');
hold on 
plot(deltaZ/lambda,RZ_NLOS,'r');
grid on
plot(deltaY/lambda,a,'m');
xlabel('distance[m]');ylabel('Correlation');
title('Spatial correlation of the new channel model for each direction in narrowband - NLOS case')
legend('X NLOS case','Y NLOS case','Z NLOS case','clarke model')
