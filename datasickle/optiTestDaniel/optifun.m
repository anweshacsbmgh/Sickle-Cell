function [ Err ] = optifun( x0,vel,tt,inten )
% load('data.mat')
[tt, ind]=sort(tt);
y = (vel(ind))';
y = smooth(y,10,'moving');
yInterp = interp1(tt(2:end),y(2:end)',1:max(tt),'spline');
x = (inten(ind))';
x = smooth(x,10,'moving');
xInterp = interp1(tt(2:end),x(2:end)',1:max(tt),'spline');
timer = tt';
xInterp = max(xInterp,0); %to remove -ve values
yInterp = max(yInterp,0); %to remove -ve values
yF = fft(yInterp);
xF = fft(xInterp);
Err=0;
tau0 = x0;
while tau0>50
    break
end

tlow = find(timer>tau0,1);
for f = 1:218 
    Err = Err+ (real(yF(f)) - real(xF(f)*exp(-1i*2*3.14*f*tau0))).^2 + (imag(yF(f)) - imag(xF(f)*exp(-1i*2*3.14*f*tau0))).^2; % Sum of square of error
end

end



