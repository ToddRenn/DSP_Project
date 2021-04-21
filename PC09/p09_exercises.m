% This file contains the lab exercises for P-09
%% 3.1: Deconvolution Experiment for 1D Filters
clc
clear
% w[n] = x[n] - 0.9x[n-1]
xx_3_1 = 256*(rem(0:100,50)<10);    % bool(#mod50 < 10) -> logical 1/0
%xx = [zeros(1) xx];
bb_3_1 = [1 -0.9];
w_n_3_1 = firfilt(bb_3_1,xx_3_1);

nn=1:75;

figure(1)
subplot(2,1,1);
stem(nn,xx_3_1(nn)); title('x[n]');
subplot(2,1,2);
stem(nn,w_n_3_1(nn),'filled'); title('w[n]');

%% 3.1.1: Restoration Filter
r = 0.9;
M = 22;

% a) Process w_n with restoration filter
yy_3_1_1 = restore(w_n_3_1,M,r);

% b) Plot the two
figure(2)
subplot(2,1,1);
stem(nn,xx_3_1(nn)); title('x[n]');
subplot(2,1,2);
stem(nn,yy_3_1_1(nn),'filled'); title('y[n]');

% c) Error estimation
error = zeros(1,length(xx_3_1));
for n = 1:length(error)
    error(n)=xx_3_1(n)-yy_3_1_1(n);
end
n=1:50;
figure(3)
plot(n,error(n));title('Error');

%% 3.1.2 Worst-Case Error
worst_case=max(error);

%% 3.1.3 An Echo Filter
% y1[n] = x1[n] + r*x1[n-P]
% a) Delay = 0.2s @ 90% given Fs = 8kHz
Fs = 8e3;
r = 0.9;  % 90%
P = Fs * 0.2; % 1600

% b) Describe filter coefficients and determine length.
% The coefficients will be 1 followed by P-2 zeros followed by 0.9 and
% the length will be P+1
bb_3_1_3 = zeros(1,P-1);
bb_3_1_3 = [1 bb_3_1_3 0.9];

% c) Implement echo filter on x2 from labdat.mat
load labdat.mat
echo=firfilt(bb_3_1_3,x2);
sound(echo,Fs);   % This is crazy distorted... but so is original x2

%% 3.2 Cascading Two Systems
% Filter 1: w[n] = x[n] - q*x[n-1]
q = 0.9;
filt1 = [1 -q];
% Filter 2: From 3.1.1 (restoration filter)
M = 22;
r = 0.9;
filt2 = zeros(1,M);
for n = 1:length(filt2)+1
    filt2(n) = 0.9^(n-1);
end
%% 3.2.1 Overall Impulse Response
% a) Implement cascade with q=0.9, r=0.9, M=22.
% Use two calls to firfilt()
y_3_2_1 = firfilt(filt2, firfilt(filt1, x1));

% Plot impulse response.
hold on
plot(x1)
plot(y_3_2_1(1:length(x1)))
title('Deconvolution');
legend('original','deconv');
hold off

% b) Work out impulse response h(n) of this cascade by hand. Verify
% correctness.

% c) What condition on h1(n) and h2(n) needs to be met for this to be
% perfect deconvolution? ANS: H1H2 = 1

%% 3.2.2 Distoring and Restoring Images
% a) Load echart.mat
load echart.mat

% b) Apply filt1 horizontally and then vertically (q=0.9)
w = firfilt(filt1, echart);    % Horizontal
ech90 = firfilt(filt1', w);    % Vertical

% c) Deconvolve with filt2 (M=22, r=0.9)
y = firfilt(filt2',ech90);      % Deconvolve vertically
y_fin = firfilt(filt2,y);   % Deconvolve horizontally
% NOTE: The order of deconvolution doesn't matter (because cascade systems
% are just multiplications in frequency domain)

figure(1); image(echart); title('Original image');
figure(2); image(ech90); title('Dist');
figure(3); image(y_fin); title('Deconv');

% Describe visual appearance of output: Thick dark line at the bottom,
% "ghost" letters appearing at the same offset below the other letters.
% Colors still inverted from original.

% Explain why ghosts exist: 

% Calculate how big the ghosts (echoes) are using worst-case error to say
% how big the ghosts are relative to "black-white" transitions which are
% 0-255: 

%% 3.2.3 A Second Restoration Experiment
% a) Deconvolve with different number of taps
% M = 11, r = 0.9
filt2_1 = fir(11,0.9);
y1 = firfilt(filt2_1',ech90);  % Deconvolve vertically
y1_fin = firfilt(filt2_1,y1);   % Deconvolve horizontally
figure(4); image(y1_fin); title('Deconv: M = 11');

% M = 22, r = 0.9
filt2_2 = fir(22,0.9);
y2 = firfilt(filt2_2',ech90);  % Deconvolve vertically
y2_fin = firfilt(filt2_2,y2);   % Deconvolve horizontally
figure(5); image(y2_fin); title('Deconv: M = 22');

% M = 33, r = 0.9
filt2_3 = fir(33,0.913);
y3 = firfilt(filt2_3',ech90);  % Deconvolve vertically
y3_fin = firfilt(filt2_3,y3);   % Deconvolve horizontally
figure(6); image(y3_fin); title('Deconv: M = 33');

% Best result was M = 33. As you increase the number of taps, the
% background color becomes more similar to original image and the ghost
% images fade into it. However, the bottom blue bar becomes larger. This is
% a result of increasing the number of rows and can be eliminated by taking
% only the first ~250 rows. 
% Explain how this is a result of the math: 
%   Hint: determine impulse response and relate it to theimage
%   Hint: use dconvdemo to generate h(n) for the systems (like the warmup)

% b) Gray scale has 256 levels. How large is the worst-case error in terms
% of gray levels? Do for each of the filters in (a).

%% Extra Credit: Filtering Music Waveform
[nggyu,Fs] = audioread('nggyu.wav');
% Range of interest: 1.915e6:2.65e6
sample = nggyu(1.915e6:2.65e6);
% Fs = 44100 Hz
a = 0.95;    % Strength of echo
T = 0.15;    % Delay (s)
P = T*Fs;   % # of samples to achieve delay

% a) Define the delay filter coefficients for 0.15s delay
delay = zeros(1,P-1);
delay = [1/(1+a) delay a/(1+a)];

% b) Filter the song file. 
%   QUESTION: Describe sound and use impulse response to
% explain why it sounds that way.
song_delay = firfilt(delay,sample);
clear sound
%sound(song_delay,Fs)    % DJ, rock dem beatz

% c) Reverberation (multiple echoes) Cascade the filters.
    % QUESTION: Derive impulse response for cascade of 4 delays
song_delay4 = firfilt(delay, firfilt(delay, firfilt(delay, firfilt(delay, sample))));
%clear sound
sound(song_delay4,Fs)

% d) QUESTION: Describe sound you hear and use impulse response to explain
% e) Plot the original + delay
figure(1); inout(sample, song_delay, 1, Fs, 4); title('Single delay');
figure(2); inout(sample, song_delay4, 1, Fs, 4); title('Reverb');

%% FUNCTIONS
function y = restore(w,M,r)
    % This is a restoration filter.
    %   w: Filtered signal
    %   M: Number of taps
    %   r: Weighted modifier
    y = zeros(1, length(w));
    for n = 1:length(w)
        for l = 0:M
            if (n-l < 1) || (n-l > length(w))
                continue
            end
            y(n) = y(n) + (r^l)*w(n-l);
        end
    end
end

function y = fir(M,r)
    % This function creates a FIR filter based on filter2 equation
    %   M: Number of taps
    %   r: Weighted modifier
    y = zeros(1,M);
    for n = 1:length(y)
        y(n) = r^(n-1);
    end
end