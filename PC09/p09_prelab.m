% This file contains the pre-lab for P-09
%% Filtering via Convolution
clc
clear

% Define b_k constants for h(n)
bb = 1/3*ones(1,3);

% Define input signal
xx = [ones(1,10),zeros(1,5)];

% Convolve
yy = firfilt(bb,xx);   % NOTE: length(y) = length(bb)+length(xx)-1

% Define discrete time points
nn = 1:length(xx);
subplot(2,1,1);
stem(nn,xx(nn))
title('x[n]')
subplot(2,1,2);
stem(nn,yy(nn),'filled')
title('y[n]')
xlabel('Time Index (n)')
