function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
global Sh
pl = 0;
ql = 1;
pr = ur;
qr = Sh;