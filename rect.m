%改进过的rect函数，可用于计算二维和一维情况下的rect

function result=rect(t)
% t=[0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9,0.8,0.7,0.6,0.5,0.3,0.2,0.1];
  [M,N]=size(t);
 length_t=M*N;
 vector=zeros(1,length_t);
  index=find(abs(t)<= 0.5);
  vector(index)=1;
  result=reshape(vector,M,N);
%   figure;plot(t);
%   figure;plot(result);
end
