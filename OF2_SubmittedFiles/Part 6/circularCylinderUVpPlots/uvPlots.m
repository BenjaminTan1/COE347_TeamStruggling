clear
close
clc

fig = 0;

U = readmatrix('U');
p = readmatrix('p');
% A = table2array(U);

fig = fig+1;
figure(fig)
hold on

plot(U(:,1),U(:,2),'blue');
plot(U(:,1),U(:,3),'red');
xlabel('t/(D/U)');
ylabel('Velocity Component (u/U, v/U)');
legend ('u/U','v/U');
title('Probe 1 Time History of u/U, v/U Sampled at (x,y)=(5.5,-0.5)');
hold off

fig = fig+1;
figure(fig)
hold on
plot(p(:,1),p(:,2),'blue');
xlabel('t/(D/U)');
ylabel('Pressure (p/rho*U^2)');
title('Probe 1 Time History of p/rho*U^2 Sampled at (x,y)=(5.5,-0.5)');

fig = fig+1;
figure(fig)
hold on
plot(U(:,1),U(:,5),'blue');
plot(U(:,1),U(:,6),'red');
xlabel('t/(D/U)');
ylabel('Velocity Component (u/U, v/U)');
legend ('u/U','v/U');
title('Probe 2 Time History of u/U, v/U Sampled at (x,y)=(5.5,+0.5)');

fig = fig+1;
figure(fig)
hold on
plot(p(:,1),p(:,3),'red');
xlabel('t/(D/U)');
ylabel('Pressure (p/rho*U^2)');
title('Probe 2 Time History of p/rho*U^2 Sampled at (x,y)=(5.5,+0.5)');

