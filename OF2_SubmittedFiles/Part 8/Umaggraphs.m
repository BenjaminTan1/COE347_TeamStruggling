clear
close
clc

fig = 0;

U = load('data.mat');

fig = fig+1;
figure(fig)
hold on
grid on
magn = sqrt((U.cc110refinedU{:,2}).^2+(U.cc110refinedU{:,3}).^2);
plot(U.cc110refinedU{:,1},magn);
xlabel('t/(D/U)');
ylabel('Velocity Magnitude');
title('Probe 1 Time History of U Sampled at (x,y)=(5.5,-0.5)');
xlim([60,130]);
hold off
