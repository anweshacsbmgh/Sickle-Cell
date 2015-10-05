clear 
close all
clc
Tau0 = 4; %initial guess for Tau
load('TS.2015.06.23.r10.2009.01.14.20150630123538.mat')

for i = 1:length(stTS)
vel(i) = stTS(i).velocity(1).*sign(stTS(i).velocity(1));
tt (i) = stTS(i).time-stTS(1).time;
inten(i) = stTS(i).intensity(1);
oo(i) = stTS(i).oxygen;
end
diff = tt(2:end)-tt(1);
diff = seconds(diff);
vel = (vel-min(vel))/(max(vel)-min(vel));
inten = (inten-min(inten))/(max(inten)-min(inten));

idlow=tt<1055;
idup=tt>1450;
vel(idup)=[];
vel(idlow)=[];
inten(idup)=[];
inten(idlow)=[];
tt(idup)=[];
tt(idlow)=[];

% options for optimization. You can change the algorithm used :
% 'interior-point', 'sqp', 'active-set', 'trust-region-reflective' when
% using fmincon
options = psoptimset('Display','iter','TolFun',1e-10,'TolX',1e-10,'TolMesh',1e-12); %options for patternsearch
% options =optimset('Display','iter','TolFun',1e-10,'TolX',1e-10,'Algorithm','sqp');

% no linear or nonlinear constraints (equality or inequality), so "[]"
% lower bound for Tau is 0 and upper bound is 25. It is actually the number
% of measurements and not explicitly the time
 

Tau = patternsearch(@(x)optifun(x,vel,tt,inten),Tau0,[],[],[],[],0,20,[],options); 
% Tau = fmincon(@optifun,Tau0,[],[],[],[],0,20,[],options);
%  Tau = fminsearch(@optifun,Tau0,options);

% % Plotting based on the optimized tau
% load('data.mat')
% [tt, ind]=sort(nlowTime);
% y = (normVel(ind))';
% y = smooth(y,0.1,'rlowess');
% yInterp = interp1(tt(2:end),y(2:end)',1:218,'spline');
% x = (normInt(ind))';
% x = smooth(x,0.1,'rlowess');
% xInterp = interp1(tt(2:end),x(2:end)',1:218,'spline');
% xInterp = max(xInterp,0); %to remove -ve values
% yInterp = max(yInterp,0); %to remove -ve values
plot(1:length(tt),vel)
hold on
% 
for i = ceil(Tau)+1:length(vel)
yfit(i-ceil(Tau)) = inten(i);
end
plot(1:length(tt)-Tau,yfit,'o');
legend('Velocity_{measured}','Velocity_{fit}')