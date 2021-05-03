%% 3 the Lab 

%% 3.1 Edge detection and Location via 1-D filters
clear 
clc
xx = 255*(rem(1:159,30)>19);

bb = [1/255, -1/255];

yy = firfilt(xx,bb);

figure
subplot(2,1,1);
stem(xx)
title('Discrete Sig')
subplot(2,1,2);
stem(yy) 
title('Edges')
xlabel("Time Index (n)")

% true value setting
d = [];

for i =1:length(yy)
    if abs(yy(i))==1
        d(i) = true;
    else 
        d(i) = false;
        
    end
end

edges = find(d)

figure
stem(edges)
title('Edges in Signal')
xlabel('nth edge')
ylabel('location of edge')

%% Bar code
clear, clc
HP = imread('HP110v3.png');

%pick a row you want out of the bar code
m = 250;
xn = HP(m,:);
%the edge detector filter
bb = [1, -1];

%do the edge detection
yn = firfilt(bb, xn);

figure
subplot(2,1,1);
stem(xn)
title('Discrete Sig')
subplot(2,1,2);
stem(yn) 
title('Edges')
xlabel("Time Index (n)")

%turn the bar code white to black shade to binary
for i =1:length(yn)
    if abs(yn(i))>=60
        d(i) = true;
    else 
        d(i) = false;
        
    end
end

%this find all the edges
L = find(d);

figure
stem(L)
title('Edges in Signal')
xlabel('nth edge')
ylabel('location of edge')


% Delta from L 
%This finds the width of each strip of color
Delta = firfilt(bb,L);
Delta(length(Delta)) = [];

%plot
figure
subplot(2,1,1);
stem(L)
title('Edge Locations')
subplot(2,1,2);
stem(Delta) 
title('Deltas')
xlabel("nth edge")

% e : logic statment

%Each bar code consists of 12 numbers
%each number is 7 units wide

%3 unit buffer on each end

%5 unit buffer in the middle

%Width = 7units*12 + 2*3units + 5 units

%Width = 95units

%% f/g/h/i for the HP110v3 Bar code

Theta = max(Delta)/4;

relativeD = round(Delta/Theta)
relativeD = relativeD(7:65);
code = decodeUPC(relativeD)



%% f/g/h/i for the OFFv3 Bar code
Theta = 11/4;

relativeD = round(Delta/Theta);
relativeD = relativeD(5:63);
code = decodeUPC(relativeD)