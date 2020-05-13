% this is the main file for ECE599 matrix analysis for SP and ML mid term
clear;
clc;
close all; 

%% signal generation

N = 15; % number of sensors
K = 5; % number of sources
T = 200; % number of samples

theta = [-60,-55,-52,13,19]*pi/180; % note that sine in matlab does not directly take spatial degrees.

A = zeros(N,K);
for n=1:N
    for k=1:K
     A(n,k)=exp(-1i*2*pi*(n-1)*(1/2)*sin(theta(k)));
    end
end


S = randn(K,T);
W = .01*randn(N,T);

X = A*S + W;


%% MUSIC

[Smusic] = MUSIC(X,K);  % the output should be the MUSIC spectrum defined over theta_scan
    
theta_scan_deg = -90:0.05:90; 
theta_scan = theta_scan_deg*pi/180;
figure(1); 
plot(theta_scan*180/pi,10*log10(Smusic),'b-','linewidth',1.5); hold on
grid on;
title('MUSIC Spectrum','fontsize',16);
xlabel('degree','fontsize',16);
ylabel('Magnitude, in dB','fontsize',16);
xlim([-90,90])

%% ESPRIT
[theta_esprit] = ESPRIT(X,K); % the output theta_esprit should be a column vectors contain all the DOAs (in degree) you estimated.

disp(['the ESPRIT estimation result are: ',num2str(theta_esprit')])

