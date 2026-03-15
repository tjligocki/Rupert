function [ x ] = normnxn(x)
nn = length(x);

n = (-1 + sqrt(1 + 8*nn)) / 2;

ii = 1;
for i = 1:n-1
  temp = x(ii:ii + n-i);
  temp = temp/norm(temp);
  x(ii:ii + n-i) = temp;
  ii = ii + n-i + 1;
end

x(nn) = x(nn-2);
