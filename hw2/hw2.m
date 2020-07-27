img1=imread('input/blurry_moon.tif');
img2=imread('input/skeleton_orig.bmp');

out1=laplacian(img1);
figure,
subplot(1,2,1);imshow(img1);title('origin images');
subplot(1,2,2);imshow(out1);title('images after laplacian');
out2=laplacian(img2);
figure,
subplot(1,2,1);imshow(img2);title('origin images');
subplot(1,2,2);imshow(out2);title('images after laplacian');

high_boost(img1);
high_boost(img2);



function out= laplacian(img)
kernel = -1 * ones(3);
kernel(2,2) = 8;  % Now kernel = [-1,-1,-1; -1,8,-1; -1,-1,-1]
sharp = imfilter(img, kernel);
out=img-sharp;
end

function []= high_boost(img)
h = fspecial('gaussian',5,2.5);
% blur the image
blurred_img = imfilter(img,h);
% subtract blurred image from original
diff_img = img - blurred_img;
% add difference to the original image
highboost_img = img + 3*diff_img;

figure,
subplot 121
imshow(img,[]);
title('Original Image')
subplot 122
imshow(highboost_img,[]);
title('HighBoosted Image')

end
