clear 
close all
Sh = 0.1;
global Sh
m = 0; % 0 for slab, 1 for cylinder, 2 for sphere
x = linspace(0,1,20);
t = linspace(0,500,1000);
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

% Extract the first solution component as u.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
surf(x,t,u) 
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance x')
ylabel('Time t')

% A solution profile can also be illuminating.
figure
plot(x,u(end,:))
title('Solution at t = t_{end}')
xlabel('Distance x')
ylabel('u(x,t_{end})')

