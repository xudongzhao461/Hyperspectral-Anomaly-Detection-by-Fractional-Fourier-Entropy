%% demo for Hyperspectral Anomaly Detection by Fractional Fourier Entropy
%--------------Brief description-------------------------------------------
%
% 
% This demo implements FRFE-RX destriping for HSI [1]
%
%
% More details in:
%
% [1] Ran Tao£¬Xudong Zhao£¬Wei Li ,Hengchao Li and Qian Du,  "Hyperspectral Anomaly Detection 
% by Fractional Fourier Entropy," IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing.
%
% contact: liwei089@ieee.org (Wei Li),zhaoxudong@bit.edu.cn (Xudong Zhao)
% 
clear all;  clc; 
close all

addpath(genpath('../Dataset'))

TIR = load('abu-airport-4.mat');
DataTest = TIR.data;
DataTest = DataTest./max(DataTest(:));
mask = double(TIR.map);
[rows cols bands] = size(DataTest);
M = rows * cols;

X = reshape(DataTest, rows*cols, bands);
A = find(mask==1);
Datatmp = reshape(DataTest, rows*cols, bands);

% Global RX 
r0 = RX(X');

% Derivative RX
X2 = abs(X(:,5:end)-X(:,1:end-4)); % step = 4
r1 = RX(X2');

% DWT-RX
for i = 1:M
    [ca,cd] = dwt(X(i, :),'db1','mode','sym');
    XW(i, :)=[ca];
end
r3 = RX(XW');

tic
disp('Running ROC...');
mask = reshape(mask, 1, M);
anomaly_map = logical(double(mask)==1);
normal_map = logical(double(mask)==0);
r_max = max(r0(:)); r_min = min(r0(:));
taus = linspace(r_min, r_max, 5000);
for index2 = 1:length(taus)
  tau = taus(index2);
  anomaly_map_rx = (r0 > tau);
  PF0(index2) = sum(anomaly_map_rx & normal_map)/sum(normal_map);
  PD0(index2) = sum(anomaly_map_rx & anomaly_map)/sum(anomaly_map);
end

r_max = max(r1(:)); r_min = min(r1(:));
taus = linspace(r_min, r_max, 5000);
for index2 = 1:length(taus)
  tau = taus(index2);
  anomaly_map_rx = (r1 > tau);
  PF1(index2) = sum(anomaly_map_rx & normal_map)/sum(normal_map);
  PD1(index2) = sum(anomaly_map_rx & anomaly_map)/sum(anomaly_map);
end

r_max = max(r3(:)); r_min = min(r3(:));
taus = linspace(r_min, r_max, 5000);
for index2 = 1:length(taus)
  tau = taus(index2);
  anomaly_map_rx = (r3 > tau);
  PF3(index2) = sum(anomaly_map_rx & normal_map)/sum(normal_map);
  PD3(index2) = sum(anomaly_map_rx & anomaly_map)/sum(anomaly_map);
end

% FrFT + RX
[FrFE,order]=FRFEorder(DataTest);
im1 = zeros(rows, cols, bands);
for i = 1:rows
    for j = 1:cols
        % amptitude
        im1(i,j,:) = Disfrft(squeeze(DataTest(i,j,:)),order);
    end
end
im1 = center_standard(abs(im1));
XT = reshape(im1, rows*cols, bands);
r2 = RX(XT');

r_max = max(r2(:)); r_min = min(r2(:));
taus = linspace(r_min, r_max, 5000);
for index2 = 1:length(taus)
  tau = taus(index2);
  anomaly_map_rx = (r2 > tau);
  PF2(index2) = sum(anomaly_map_rx & normal_map)/sum(normal_map);
  PD2(index2) = sum(anomaly_map_rx & anomaly_map)/sum(anomaly_map);
end


% DFT + RX
im1 = zeros(rows, cols, bands);
for i = 1:rows
    for j = 1:cols
        % amptitude
        im1(i,j,:) = fft(squeeze(DataTest(i,j,:)));
    end
end
im1 = center_standard(abs(im1));
XT = reshape(im1, rows*cols, bands);
r5 = RX(XT');

r_max = max(r5(:)); r_min = min(r5(:));
taus = linspace(r_min, r_max, 5000);
for index2 = 1:length(taus)
  tau = taus(index2);
  anomaly_map_rx = (r5 > tau);
  PF5(index2) = sum(anomaly_map_rx & normal_map)/sum(normal_map);
  PD5(index2) = sum(anomaly_map_rx & anomaly_map)/sum(anomaly_map);
end


area0 = sum((PF0(1:end-1)-PF0(2:end)).*(PD0(2:end)+PD0(1:end-1))/2);
area1 = sum((PF1(1:end-1)-PF1(2:end)).*(PD1(2:end)+PD1(1:end-1))/2);
area3 = sum((PF3(1:end-1)-PF3(2:end)).*(PD3(2:end)+PD3(1:end-1))/2);
area2 = sum((PF2(1:end-1)-PF2(2:end)).*(PD2(2:end)+PD2(1:end-1))/2);
area5 = sum((PF5(1:end-1)-PF5(2:end)).*(PD5(2:end)+PD5(1:end-1))/2);

figure,
plot(PF0, PD0, 'r-.', 'LineWidth', 2);  hold on
plot(PF5, PD5, 'm-', 'LineWidth', 2);  hold on
plot(PF3, PD3, 'c-', 'LineWidth', 2);  hold on
plot(PF1, PD1, 'b--', 'LineWidth', 2);  hold on
plot(PF2, PD2, 'k-', 'LineWidth', 2);  grid on
xlabel('Probability of false alarm'); ylabel('Probability of detection');
legend('Global-RX','DFT-RX','DWT-RX', 'Deriv-RX','FrFT-RX')
axis([0 0.3 0.8 1])  ;