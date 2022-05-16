%% Read Files
clc;
% 25
filename = "25\25_pTop.xlsx";
pTop25 = readtable(filename);
pTop25 = table2array(pTop25);
filename = "25\25_pBot.xlsx";
pBot25 = readtable(filename);
pBot25 = table2array(pBot25);
filename = "25\25_pOutlet.xlsx";
pOut25 = readtable(filename);
pOut25 = table2array(pOut25);

% 50
filename = "50\50_pTop.xlsx";
pTop50 = readtable(filename);
pTop50 = table2array(pTop50);
filename = "50\50_pBot.xlsx";
pBot50 = readtable(filename);
pBot50 = table2array(pBot50);
filename = "50\50_pOutlet.xlsx";
pOut50 = readtable(filename);
pOut50 = table2array(pOut50);

% 100
filename = "100\100_pTop.xlsx";
pTop100 = readtable(filename);
pTop100 = table2array(pTop100);
filename = "100\100_pBot.xlsx";
pBot100 = readtable(filename);
pBot100 = table2array(pBot100);
filename = "100\100_pOutlet.xlsx";
pOut100 = readtable(filename);
pOut100 = table2array(pOut100);

% Box
filename = "Box\Box_pTop.xlsx";
pTopBox = readtable(filename);
pTopBox = table2array(pTopBox);
filename = "Box\Box_pBot.xlsx";
pBotBox = readtable(filename);
pBotBox = table2array(pBotBox);
filename = "Box\Box_pOutlet.xlsx";
pOutBox = readtable(filename);
pOutBox = table2array(pOutBox);

%% Figures
% Top
figure;
plot(pTop25(:,1),pTop25(:,4));
hold on;
plot(pTop50(:,1),pTop50(:,4));
hold on;
plot(pTop100(:,1),pTop100(:,4));
hold on;
scatter(pTopBox(:,1),pTopBox(:,4));
xlim([-400,400]);
legend('Coarse Spike','Medium Spike','Refined Spike', 'Refined Box');
title('Pressure at Top');
xlabel('x (m)')
ylabel('pressure (Pa)')

% Bot
figure;
plot(pBot25(:,1),pBot25(:,4));
hold on;
plot(pBot50(:,1),pBot50(:,4));
hold on;
plot(pBot100(:,1),pBot100(:,4));
hold on;
scatter(pBotBox(:,1),pBotBox(:,4));
xlim([-400,400]);
legend('Coarse Spike','Medium Spike','Refined Spike', 'Refined Box');
title('Pressure at Bottom');
xlabel('x (m)')
ylabel('pressure (Pa)')

% Outlet
figure;
plot(pOut25(:,2),pOut25(:,4));
hold on;
plot(pOut50(:,2),pOut50(:,4));
hold on;
plot(pOut100(:,2),pOut100(:,4));
hold on;
scatter(pOutBox(:,2),pOutBox(:,4));
legend('Coarse Spike','Medium Spike','Refined Spike', 'Refined Box');
xlim([-60,60]);
title('Pressure in Wake');
xlabel('y (m)')
ylabel('pressure (Pa)')

