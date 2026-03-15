function [ x,y,z,v ] = getdata3( file )
%GETDATA3 Read in 3D angular data
%  Read in a file which is N^3 by 4 and extract the columns (z,y,x,v) and
%  reshape them as NxNxN matrices
alldata = load(file);
num = size(alldata);
n = round(num(1)^(1.0/3.0));
z = reshape(alldata(:,1),n,n,n);
y = reshape(alldata(:,2),n,n,n);
x = reshape(alldata(:,3),n,n,n);
v = reshape(alldata(:,4),n,n,n);
min(min(min(v)))
max(max(max(v)))