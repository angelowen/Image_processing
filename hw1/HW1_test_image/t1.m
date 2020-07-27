
split16('Lena.bmp',1);
split16('Lena.bmp',2);


function []=split16(filename,op)
A=imread(filename);
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

end

