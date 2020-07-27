
clear all;
close all;
img1 = imread('images/pool.png');
img2 = imread('images/peppers.png');
img3 = imread('images/baboon.png');
mysobel(img1,30,1)
mysobel(img2,50,2)
mysobel(img3,100,3)

function []=mysobel(img,threshold,idx)
figure(idx);
subplot(1,2,1);imshow(img);title('input image');
img = rgb2gray(img);
C = double(img);
for i = 1:size(C,1)-2
    for j = 1:size(C,2)-2
        x = ((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        y = ((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
        
        img(i,j) = sqrt(x.^2+y.^2);
        
    end
end
img = max(img,threshold);
img(img==round(threshold))=0;
img = uint8(img);
name=['Edges Detected after threshold ' int2str(threshold)];
subplot(1,2,2);imshow(~img);title(name);
end