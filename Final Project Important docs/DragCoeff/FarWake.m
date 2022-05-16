filename = "Spike\neg1500_Spike_p.xlsx";
pos1500_Spike_p = readtable(filename);
pos1500_Spike_p = table2array(pos1500_Spike_p);
filename = "Box\neg1500_Box_p.xlsx";
pos1500_Box_p = readtable(filename);
pos1500_Box_p = table2array(pos1500_Box_p);

figure;
plot(pos1500_Spike_p(:,2),pos1500_Spike_p(:,4));
hold on;
plot(pos1500_Box_p(:,2),pos1500_Box_p(:,4));
legend('Refined Spike', 'Refined Box');
title('Pressure at Inlet');
xlabel('y (m)')
ylabel('pressure (Pa)')