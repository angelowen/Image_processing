%% restart
close all;
clc;
A=imread('Lena.bmp');
Img=A;
%% WINDOW SIZE
M=16;
N=16;
mid_val=round((M*N)/2);
%% FIND THE NUMBER OF ROWS AND COLUMNS TO BE PADDED WITH ZERO
in=0;
for i=1:M
    for j=1:N
        in=in+1;
        if(in==mid_val)
            PadM=i-1;
            PadN=j-1;
            break;
        end
    end
end
%% PADDING THE IMAGE WITH ZERO ON ALL SIDES
B=padarray(A,[PadM,PadN]);

for i= 1:size(B,1)-((PadM*2)+1)
    for j=1:size(B,2)-((PadN*2)+1)
      
        cdf=zeros(256,1);
        inc=1;
        for x=1:M
            for y=1:N
                
  %FIND THE MIDDLE ELEMENT IN THE WINDOW          
                if(inc==mid_val)
                    ele=B(i+x-1,j+y-1)+1;
                end
                    pos=B(i+x-1,j+y-1)+1;
                    cdf(pos)=cdf(pos)+1;
                   inc=inc+1;
            end
        end  
        %COMPUTE THE CDF FOR THE VALUES IN THE WINDOW
        for l=2:256
            cdf(l)=cdf(l)+cdf(l-1);
        end
            Img(i,j)=round(cdf(ele)/(M*N)*255); 
         
    end
end
%% draw img
figure,
subplot(2,1,1);imshow(A);title('origin images');
subplot(2,1,2);imshow(Img);title('images after Histogram Equalization');

split16(A,1,'origin images split 16');
split16(A,2,'Histogram of origin images split 16');
split16(Img,1,'16 split images after Histogram Equalization');
split16(Img,2,' Histogram of 16 split images after Histogram Equalization');

figure,
subplot(2,1,1);histogram(A);title('Before Local Histogram Equalization');
subplot(2,1,2);histogram(Img);title('After Local Histogram Equalization'); 
%% Second image
A=imread('Peppers.bmp');
Img=A;
B=padarray(A,[PadM,PadN]);

for i= 1:size(B,1)-((PadM*2)+1)
    for j=1:size(B,2)-((PadN*2)+1)
      
        cdf=zeros(256,1);
        inc=1;
        for x=1:M
            for y=1:N
                
  %FIND THE MIDDLE ELEMENT IN THE WINDOW          
                if(inc==mid_val)
                    ele=B(i+x-1,j+y-1)+1;
                end
                    pos=B(i+x-1,j+y-1)+1;
                    cdf(pos)=cdf(pos)+1;
                   inc=inc+1;
            end
        end  
        %COMPUTE THE CDF FOR THE VALUES IN THE WINDOW
        for l=2:256
            cdf(l)=cdf(l)+cdf(l-1);
        end
            Img(i,j)=round(cdf(ele)/(M*N)*255); 
         
    end
end
% draw img
figure,
subplot(2,1,1);imshow(A);title('origin images');
subplot(2,1,2);imshow(Img);title('images after Histogram Equalization');
split16(A,1,'origin images split 16');
split16(A,2,'Histogram of origin images split 16');
split16(Img,1,'16 split images after Histogram Equalization');
split16(Img,2,' Histogram of 16 split images after Histogram Equalization');

figure,
subplot(2,1,1);histogram(A);title('Before Local Histogram Equalization');
subplot(2,1,2);histogram(Img);title('After Local Histogram Equalization'); 

%% split16

function []=split16(A,op,titlename)

[m,n]=size(A);
m1=fix(m/4)-1;
n1=fix(n/4)-1;

id=1;
figure() ;
for i=1:m1:m-m1
    for j=1:n1:n-n1
        subplot(4, 4, id) ;
        if op==1
            imshow(A(i+1:i+m1,j:j+n1,:)); 
        elseif op==2
            histogram(A(i+1:i+m1,j:j+n1,:));
        end
        id=id+1;
    end
end
sgtitle(titlename);
end













