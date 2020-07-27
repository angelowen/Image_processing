
%% start
clear all;
close all;
A = imread('images/aloe.jpg');
A1 = imread('images/church.jpg');
A2 = imread('images/house.jpg');
A3 = imread('images/kitchen.jpg');
enhancement(A,1);
enhancement(A1,2);
enhancement(A2,3);
enhancement(A3,4);
%% main
function []= enhancement(pic,num)
    Adoub = im2double(pic);
    %RGB
    [R,G,B] = imsplit(Adoub);
    Rsharpened = myHistogram(R);
    Gsharpened = myHistogram(G);
    Bsharpened = myHistogram(B);
    rgbImage = cat(3, Rsharpened,Gsharpened,Bsharpened);

    %HSI
    [H,S,I]= myrgb2hsi(pic);
    I = myHistogram(I);%histeq
    S = myHistogram(S);
    P = myhsi2rgb(H,S,I);

    %L*a*b*
    A2 = myrgb2lab(pic);
    max_luminosity = 100;
    L = A2(:,:,1)/max_luminosity;
    po = A2;
    po(:,:,1) = myHistogram(L)*max_luminosity;
    po = mylab2rgb(po);
    
    figure(num),
    subplot(2,2,1);imshow(pic); title('Original image');
    subplot(2,2,2);imshow(rgbImage); title('RGB hist');
    subplot(2,2,3);imshow(P); title('HSI result');
    subplot(2,2,4);imshow(po); title('L*a*b* histeq');
end
%% lab2rgb
function [R, G, B] = mylab2rgb(L, a, b)
if (nargin == 1)
  b = L(:,:,3);
  a = L(:,:,2);
  L = L(:,:,1);
end
% Thresholds
T1 = 0.008856;
T2 = 0.206893;
[M, N] = size(L);
s = M * N;
L = reshape(L, 1, s);
a = reshape(a, 1, s);
b = reshape(b, 1, s);
% Compute Y
fY = ((L + 16) / 116) .^ 3;
YT = fY > T1;
fY = (~YT) .* (L / 903.3) + YT .* fY;
Y = fY;
% Alter fY slightly for further calculations
fY = YT .* (fY .^ (1/3)) + (~YT) .* (7.787 .* fY + 16/116);
% Compute X
fX = a / 500 + fY;
XT = fX > T2;
X = (XT .* (fX .^ 3) + (~XT) .* ((fX - 16/116) / 7.787));
% Compute Z
fZ = fY - b / 200;
ZT = fZ > T2;
Z = (ZT .* (fZ .^ 3) + (~ZT) .* ((fZ - 16/116) / 7.787));
X = X * 0.950456;
Z = Z * 1.088754;
MAT = [ 3.240479 -1.537150 -0.498535;
       -0.969256  1.875992  0.041556;
        0.055648 -0.204043  1.057311];
RGB = max(min(MAT * [X; Y; Z], 1), 0);
R = reshape(RGB(1,:), M, N) * 255;
G = reshape(RGB(2,:), M, N) * 255;
B = reshape(RGB(3,:), M, N) * 255; 
if ((nargout == 1) | (nargout == 0))
  R = uint8(round(cat(3,R,G,B)));
end
end
%% rgb2lab
function [L,a,b] = myrgb2lab(R,G,B)
if nargin == 1
  B = double(R(:,:,3));
  G = double(R(:,:,2));
  R = double(R(:,:,1));
end
if max(max(R)) > 1.0 || max(max(G)) > 1.0 || max(max(B)) > 1.0
  R = double(R) / 255;
  G = double(G) / 255;
  B = double(B) / 255;
end
% Set a threshold
T = 0.008856;
[M, N] = size(R);
s = M * N;
RGB = [reshape(R,1,s); reshape(G,1,s); reshape(B,1,s)];
% RGB to XYZ
MAT = [0.412453 0.357580 0.180423;
       0.212671 0.715160 0.072169;
       0.019334 0.119193 0.950227];
XYZ = MAT * RGB;
% Normalize for D65 white point
X = XYZ(1,:) / 0.950456;
Y = XYZ(2,:);
Z = XYZ(3,:) / 1.088754;
XT = X > T;
YT = Y > T;
ZT = Z > T;
Y3 = Y.^(1/3); 
fX = XT .* X.^(1/3) + (~XT) .* (7.787 .* X + 16/116);
fY = YT .* Y3 + (~YT) .* (7.787 .* Y + 16/116);
fZ = ZT .* Z.^(1/3) + (~ZT) .* (7.787 .* Z + 16/116);
L = reshape(YT .* (116 * Y3 - 16.0) + (~YT) .* (903.3 * Y), M, N);
a = reshape(500 * (fX - fY), M, N);
b = reshape(200 * (fY - fZ), M, N);
if nargout < 2
  L = cat(3,L,a,b);
end
end
%% hsi to rgb
function rgb = myhsi2rgb(H,S,I)
    hsi = cat(3, H, S, I);
    % Extract the individual HSI component images. 
    H=hsi(:,:,1)*2*pi;
    S=hsi(:,:,2);
    I=hsi(:,:,3);
    R=zeros(size(hsi,1),size(hsi,2));
    G=zeros(size(hsi,1),size(hsi,2));
    B=zeros(size(hsi,1),size(hsi,2));
    idx=find((0<=H)&(H<2*pi/3));
    B(idx)=I(idx).*(1-S(idx));
    R(idx)=I(idx).*(1+S(idx).*cos(H(idx))./cos(pi/3-H(idx)));
    G(idx)=3*I(idx)-(R(idx)+B(idx));
    idx=find((2*pi/3<=H)&(H<4*pi/3));
    R(idx)=I(idx).*(1-S(idx));
    G(idx)=I(idx).*(1+S(idx).*cos(H(idx)-2*pi/3)./cos(pi-H(idx)));
    B(idx)=3*I(idx)-(R(idx)+G(idx));
    idx=find((4*pi/3<=H)&(H<=2*pi));
    G(idx)=I(idx).*(1-S(idx));
    B(idx)=I(idx).*(1+S(idx).*cos(H(idx)-4*pi/3)./cos(5*pi/3-H(idx)));
    R(idx)=3*I(idx)-(G(idx)+B(idx));
    rgb=cat(3,R,G,B);
    rgb=max(min(rgb,1),0);
end

%% histogram
function [res] = myHistogram(channel)

Hist=[];
for i = 1:10
    for j = 1:10
        %   Calculate the total number of excess pixels.
        NrExcess = 0;
        for nr = 1:128
            Hist(i,j,nr)=255;
            excess=Hist(i,j,nr) - 128;
            if excess > 0
                NrExcess = NrExcess + excess;
            end
        end
        %  Clip histogram and redistribute excess pixels in each bin
        binIncr = NrExcess / 128;
        upper = 128 - binIncr;
        for nr = 1:10
            if Hist(i,j,nr) > 128
                Hist(i,j,nr) = 128;
            else
                if Hist(i,j,nr) > upper
                    NrExcess = NrExcess + upper - Hist(i,j,nr);
                    Hist(i,j,nr) = 128;
                else
                    NrExcess = NrExcess - binIncr;
                    Hist(i,j,nr) = Hist(i,j,nr) + binIncr;
                end
            end
        end
        
        if NrExcess > 0
            stepSize = max(1,fix(1+NrExcess/128));
            for nr = 1:128
                NrExcess = NrExcess - stepSize;
                Hist(i,j,nr) = Hist(i,j,nr) + stepSize;
                if NrExcess < 1
                    break;
                end
            end
        end
        
    end
end
    res=adapthisteq(channel);
end


%% rgb to hsi
function [H,S,I] = myrgb2hsi(rgb) 
    rgb = im2double(rgb);  
    r = rgb(:, :, 1);  
    g = rgb(:, :, 2);  
    b = rgb(:, :, 3);    
    num = 0.5*((r - g) + (r - b));  
    den = sqrt((r - g).^2 + (r - b).*(g - b));  
    theta = acos(num./(den + eps)); %防止除數=0   
    H = theta;  
    H(b > g) = 2*pi - H(b > g);  
    H = H/(2*pi);    
    num = min(min(r, g), b);  
    den = r + g + b;  
    den(den == 0) = eps; %防止除數=0  
    S = 1 - 3.* num./den;  
    H(S == 0) = 0;  
    I = (r + g + b)/3; 
    % Combine all three results into an hsi image. 
end

