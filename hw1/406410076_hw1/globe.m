%% restart
close all;
clc;
%% read Image and Display
pic1 = imread('Lena.bmp'); %read image
%% Histogram equalization without matlab function
after_img1 = my_hist(pic1);

%% read Image and Display  img2
pic2 = imread('Peppers.bmp'); %read image
%% Histogram equalization without matlab function
after_img2 = my_hist(pic2);

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

function output = my_hist(I)
    [r, c] = size(I);
    n = r * c;
   
    output = uint8(zeros(r, c));
    f = zeros(256, 1);
    pdf = zeros(256, 1);
    cdf = zeros(256, 1);
    out = zeros(256, 1);
    cum = zeros(256, 1);

    for i = 1:r
        for j = 1:c
            value = I(i, j);
            f(value + 1) = f(value + 1) + 1;
            pdf(value + 1) = f(value + 1) / n;
        end
    end

    sum = 0;
    L = 255;
    for i = 1:size(pdf)
        sum = sum + f(i);
        cum(i) = sum;
        cdf(i) = cum(i) / n;
        out(i) = round(cdf(i) * L);
    end

    for i = 1:r
        for j = 1:c
            output(i, j) = out(I(i, j) + 1);
        end
    end
end


