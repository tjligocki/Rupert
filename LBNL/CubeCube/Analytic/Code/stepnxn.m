function [ y,err,dx ] = stepnxn(x)
toler = 1e-12;
[a,b] = trynxn(x);
dx = (a\b')';
err = norm(dx);
if err > norm(x)
  dx = dx/err;
end
y = x + dx;
y = normnxn(y);
