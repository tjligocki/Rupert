function [ a,b ] = trynxn(x)
nn = length(x);

n = (-1 + sqrt(1 + 8*nn)) / 2;

a = zeros(nn);

r = 1;
c = 1;
for i = 1:n-1
  for j = i:n
    a(r,c) = 2*x(c);
    c = c+1;
  end
  r = r+1;
end

a(r,1) = 2*x(1);

for i = 2:n
  a(r,i) = -2*x(i);
end

r = r+1;

ii = 1;
iii = n;
for i = 1:n-2
  jj = ii + iii;
  jjj = iii - 1;
  for j = i+1:n-1
    for k = j:n
      a(r,ii+k-i) = x(jj+k-j);
      a(r,jj+k-j) = x(ii+k-i);
    end
    
    jj = jj + jjj;
    jjj = jjj - 1;
    
    r = r+1;
  end
  ii = ii + iii;
  iii = iii - 1;
end

ii = 1;
iii = n;
for i = 1:n-1
  abssum = 0;
  for j = i:n-1
    abssum = abssum + abs(x(ii+j-i));
  end
  for j = i:n-1
    a(r,ii+j-i) = 2 * sign(x(ii+j-i)) * abssum;
  end
  a(r,nn) = - 2 * sign(x(nn)) * abs(x(nn));
  if a(r,nn) == 0
    a(r,nn) = 1;
  end

  ii = ii + iii;
  iii = iii - 1;
  
  r = r+1;
end

r = 1;
ii = 1;

for i=1:n-1
  b(r) = 1;
  for j=i:n
    b(r) = b(r) - x(ii)*x(ii);
    ii = ii+1;
  end
  r = r+1;
end

b(r) = -x(1)*x(1);

for i = 2:n
  b(r) = b(r) + x(i)*x(i);
end

r = r+1;

ii = 1;
iii = n;
for i = 1:n-2
  jj = ii + iii;
  jjj = iii - 1;
  for j = i+1:n-1
    b(r) = 0;
    for k = j:n
      b(r) = b(r) - x(ii+k-i)*x(jj+k-j);
    end
    
    jj = jj + jjj;
    jjj = jjj - 1;
    
    r = r+1;
  end
  ii = ii + iii;
  iii = iii - 1;
end

ii = 1;
iii = n;
for i = 1:n-1
  b(r) = 0;
  for j = i:n-1
    b(r) = b(r) + abs(x(ii+j-i));
  end
  b(r) = -b(r)^2 + abs(x(nn))^2;
  
  ii = ii + iii;
  iii = iii - 1;
  
  r = r+1;
end
