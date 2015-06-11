clear 
close all
% global Bi
m = 0; % 0 for slab, 1 for cylinder, 2 for sphere
cinf = 0.04;
x = linspace(0,1,20);
t = linspace(0,20,200);
c = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
conc(1,:) = 0.2*ones(length(c(1,:)),1);
for i = 2:length(c(:,1))
    conc(i,:) = (conc(i-1,end)-cinf)*c(i,end)+cinf;
end
% A surface plot is often a good way to study a solution.
surf(x,t,conc) 
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance x')
ylabel('Time t')
colorbar
% % A solution profile can also be illuminating.
% figure
% plot(x,conc(end,:))
% title('Solution at t = t_{end}')
% xlabel('Distance x')
% ylabel('c(x,t_{end})')

