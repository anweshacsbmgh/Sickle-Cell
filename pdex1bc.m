function [pl,ql,pr,qr] = pdex1bc(xl,cl,xr,cr,t)
Bi = 1e-2;
pl = 0;
ql = 1;
pr = Bi*cr;
qr = 1;
