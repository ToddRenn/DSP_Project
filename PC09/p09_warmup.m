% This file contains the warm-up for P-09
%% 2.1: This section requires non-Linux based GUI (con2dis)
%% 2.2: This section requires non-Linux based GUI (dconvdemo)
%% 2.3 Loading Data
%   x1: Stair-step input signal
%   x2: Speech waveform ("Oak is strong") - f = 8k samples/sec
%   xtv: Actual scan line from a digital image
%   h1: FIR coefficients
%   h2: FIR coefficients (more)
clc
clear
load labdat.mat

%% 2.4 Filtering a Signal
% a) 5-point averager on x1
bb = 1/5*ones(1,5);
y1 = firfilt(bb,x1);

% b) Plot x1, y1
nn = 1:length(x1);
subplot(2,1,1)
stem(nn,x1(nn));
title('x[n]');xlabel('n');
subplot(2,1,2);
stem(nn,y1(nn),'filled');
title('y[n]');xlabel('n');

% c) Limited x-axis plot
xmid = length(x1)/2
nn = xmid-30:xmid+30;   % 30 points from middle of signal
subplot(2,1,1)
stem(nn,x1(nn));
title('x[n]');xlabel('n');
subplot(2,1,2);
stem(nn,y1(nn),'filled');
title('y[n]');xlabel('n');

% d) Explain effect of smoother

%% 2.5 Filtering Images: 2D Convolution using conv2()
clc
clear
% a) First-difference image filtering: HORIZONTAL DIRECTION
% NOTE: this will produce an color-inverted image
load echart.mat
bdiffh = [1 -1];    % first-difference filter
yy1 = conv2(echart, bdiffh);

figure(1); image(echart); title('Original image');
figure(2); image(yy1); title('First-Difference (Hoirzontal) filtered image');

% b) First-difference image filtering: VERTICAL DIRECTION
yy2 = conv2(echart, bdiffh');

figure(3); image(yy2); title('First-Difference (Vertical) filtered image');