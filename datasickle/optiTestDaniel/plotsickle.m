clear 
close all
clc
load('TS.2015.06.23.r10.2009.01.14.20150630123538.mat')
% stTS = stTimeSeries;
for i = 1:length(stTS)
vel(i) = stTS(i).velocity(1).*sign(stTS(i).velocity(1));
tt (i) = stTS(i).time;
inten(i) = stTS(i).intensity(1);
oo(i) = stTS(i).oxygen;
end
diff = tt(2:end)-tt(1);
diff = seconds(diff);
plot(tt,oo)
figure
plot(tt,vel,'+')
figure
plot(diff,oo(2:end))
figure
plot(diff,vel(2:end),'+')
figure
plot(diff,inten(2:end),'.')
