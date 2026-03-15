function [y,errhist] = runnxn(x,cmax,cmax2)
toler = 1e-12;
err = 1;
y = x;
count = 1;
count2 = 1;
gtotal = 0;
print = 50;
emin = 10;
errhist2 = [];
while (err > toler) & (count2 <= cmax2)
  errhist = [];
  while (err > toler) & (count <= cmax)
    [y,err,dx] = stepnxn(y);
    errhist = [errhist err];
    emin = min(err,emin);
    count = count + 1;
    gtotal = gtotal + 1;
  end
  errhist2 = [errhist2 errhist];
  if mod(count2,print) == 0
    disp(sprintf('%5d: %16.10e  %16.10e',count2,min(errhist2),emin));
    errhist2 = [];
  end
  if err > toler
    pert = 0.1*(2*rand(1,length(y)) - 1);
    y = y + pert;
    y = normnxn(y);
    count = 1;
  end
  count2 = count2 + 1;
end
for i=1:5
  [y,err,dx] = stepnxn(y);
  errhist = [errhist err];
  emin = min(err,emin);
  count = count + 1;
  gtotal = gtotal + 1;
end
errhist2 = [errhist2 errhist];
disp(sprintf('%5d: %16.10e  %16.10e  %3d  %7d',count2,min(errhist2),emin,count-1,gtotal));
