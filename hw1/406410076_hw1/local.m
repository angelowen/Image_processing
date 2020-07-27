%% restart
close all;
clc;
A=imread('Lena.bmp');
B=imread('Peppers.bmp');

figure,
subplot(2,1,1);imshow(A);title('origin images');
subplot(2,1,2);imshow('Lena-local.bmp');title('images after Histogram Equalization');

split16(A,1)
split16(A,2)
split16(A,3)
split16(A,4)
figure,
subplot(2,1,1);imshow(B);title('origin images');
subplot(2,1,2);imshow('Peppers-local.bmp');title('images after Histogram Equalization');

split16(B,1)
split16(B,2)
split16(B,3)
split16(B,4)

        
function []=split16(A,op)

[m,n]=size(A);
m1=fix(m/4)-1;
n1=fix(n/4)-1;

id=1;
figure() ;
for i=1:m1:m-m1
    for j=1:n1:n-n1
        subplot(4, 4, id) ;
        if     op==1
            imshow(A(i+1:i+m1,j:j+n1,:));
        elseif op==2
            histogram(A(i+1:i+m1,j:j+n1,:));
        elseif op==3
            imshow(my_hist(A(i+1:i+m1,j:j+n1,:))); 
        elseif op==4
            histogram(my_hist(A(i+1:i+m1,j:j+n1,:)));
        
        end
        id=id+1;
    end
end

end

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