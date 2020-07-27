%% restart
close all;
clc;
%% read Image and Display
pic1 = imread('Lena.bmp'); %read image
%% Histogram equalization without matlab function
I = pic1; 

% size of the input image.
[r,c] = size(I); 

% uint 8 blamk canvas
after_img1 = uint8(zeros(r,c));

% number of pixels.
n = r*c;
% frequency pdf and cdf and some variables
f = zeros(256,1);
pdf = zeros(256,1);
cdf = zeros(256,1);
out = zeros(256,1);
cum = zeros(256,1);
% Nested for loop
for i = 1:r
    for j = 1:c
        value = I(i,j);
        f(value+1) = f(value+1)+1;
        pdf(value+1) = f(value+1)/n;
    end
end

% finding cdf
sum = 0;
L = 255;

for i = 1:size(pdf);
    sum = sum + f(i);
    cum(i) = sum;
    cdf(i) = cum(i)/n;
    out(i) = round(cdf(i)*L);
end

for i = 1:r;
    for j = 1:c;
        after_img1(i,j) = out(I(i,j)+1);
    end
end
%% img2

%% read Image and Display
pic2 = imread('Peppers.bmp'); %read image
%% Histogram equalization without matlab function
I = pic2; 

% size of the input image.
[r,c] = size(I); 

% uint 8 blamk canvas
after_img2 = uint8(zeros(r,c));

% number of pixels.
n = r*c;
% frequency pdf and cdf and some variables
f = zeros(256,1);
pdf = zeros(256,1);
cdf = zeros(256,1);
out = zeros(256,1);
cum = zeros(256,1);
% Nested for loop
for i = 1:r
    for j = 1:c
        value = I(i,j);
        f(value+1) = f(value+1)+1;
        pdf(value+1) = f(value+1)/n;
    end
end

% finding cdf
sum = 0;
L = 255;

for i = 1:size(pdf);
    sum = sum + f(i);
    cum(i) = sum;
    cdf(i) = cum(i)/n;
    out(i) = round(cdf(i)*L);
end

for i = 1:r;
    for j = 1:c;
        after_img2(i,j) = out(I(i,j)+1);
    end
end


%% draw image
figure,
subplot(2,2,1);imshow(pic1);title('origin images - Lena');
subplot(2,2,2);imshow(after_img1);title('My hist Image - Lena');
subplot(2,2,3);imshow(pic2);title('origin images - Peppers');
subplot(2,2,4);imshow(after_img2);title('My hist Image - Peppers');
figure,
subplot(2,2,1); histogram(pic1);title('Before  - Lena');
subplot(2,2,2);histogram(after_img1);title('My funtion Histogram - Lena'); 
subplot(2,2,3); histogram(pic2);title('Before  - Peppers');
subplot(2,2,4);histogram(after_img2);title('My funtion Histogram - Peppers'); 



