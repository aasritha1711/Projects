clear;
clc;
c = imread('image.jpg');
                  message = 'abcdefgh';
                  message = strtrim(message);
                  m = length(message) * 8;
                  AsciiCode = uint8(message); 
                  binaryString = transpose(dec2bin(AsciiCode,8));
                  binaryString = binaryString(:);
                  N = length(binaryString);
                  b = zeros(N,1); %b is a vector of bits
                  
%Image Encryption
                  for k = 1:N
                      if(binaryString(k) == '1')
                          b(k) = 1;
                      else
                          b(k) = 0;
                      end
                  end
                  
                  s = c;
                  height = size(c,1);
                  width = size(c,2);
                  k = 1; Array=[];l=1;my=1;
                  
                  for i = 1 : height
                      for j = 1 : width
                           LSB = mod(double(c(i,j)), 2);
                          if (k>m || LSB == b(k))
                              s(i,j) = c(i,j);
                              l=k+1;
                          else 
                              if(LSB == 1)
                                  s(i,j) = c(i,j) - 1;
                              else
                                  s(i,j) = c(i,j) + 1;
                              end
                              
                          Array(my)=l;
                          l=l+1;
                          my= my + 1;
                          end
                          k=k+1;
                      end
                  end
                   imwrite(s, 'hiddenmsgimage.jpg');
                   
                   
%Retrieval of Secret Message
 k = 1;
for i = 1 : height
  for j = 1 : width
    if( k<=m )
         b(k) = mod(double(s(i,j)),2);
    end
     k = k + 1;
  end
end

%Retrieval of Original Image
for my = 1 : numel(Array)
      temp=Array(my);
      q=(temp/width)+1;
      q=uint8(q);
      r=mod(temp,width);
      if(mod(double(s(q,r)),2)==0)
          s(q,r)=s(q,r)-1;
      else
          s(q,r)=s(q,r)+1;
      end
end

imwrite(s,'retrieved.jpg');

v=reshape(b,8,length(message));

size(b);
bin=[128 64 32 16 8 4 2 1];
v=(bin*v);
display('The Secret Message is:');
display(char(v));
