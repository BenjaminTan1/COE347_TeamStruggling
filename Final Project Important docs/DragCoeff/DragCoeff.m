%% Raw Data
clc;
clear all;
% Spikes
% Positive
filename = "Spike\pos1500_Spike_U.xlsx";
pos1500_Spike_U = readtable(filename);
pos1500_Spike_U = table2array(pos1500_Spike_U);
filename = "Spike\pos1499_Spike_U.xlsx";
pos1499_Spike_U = readtable(filename);
pos1499_Spike_U = table2array(pos1499_Spike_U);
filename = "Spike\pos1500_Spike_p.xlsx";
pos1500_Spike_p = readtable(filename);
pos1500_Spike_p = table2array(pos1500_Spike_p);

% Negative
filename = "Spike\neg1500_Spike_U.xlsx";
neg1500_Spike_U = readtable(filename);
neg1500_Spike_U = table2array(neg1500_Spike_U);
filename = "Spike\neg1499_Spike_U.xlsx";
neg1499_Spike_U = readtable(filename);
neg1499_Spike_U = table2array(neg1499_Spike_U);
filename = "Spike\neg1500_Spike_p.xlsx";
neg1500_Spike_p = readtable(filename);
neg1500_Spike_p = table2array(neg1500_Spike_p);

% Box
% Positive
filename = "Box\pos1500_Box_U.xlsx";
pos1500_Box_U = readtable(filename);
pos1500_Box_U = table2array(pos1500_Box_U);
filename = "Box\pos1499_Box_U.xlsx";
pos1499_Box_U = readtable(filename);
pos1499_Box_U = table2array(pos1499_Box_U);
filename = "Box\pos1500_Box_p.xlsx";
pos1500_Box_p = readtable(filename);
pos1500_Box_p = table2array(pos1500_Box_p);

% Negative
filename = "Box\neg1500_Box_U.xlsx";
neg1500_Box_U = readtable(filename);
neg1500_Box_U = table2array(neg1500_Box_U);
filename = "Box\neg1499_Box_U.xlsx";
neg1499_Box_U = readtable(filename);
neg1499_Box_U = table2array(neg1499_Box_U);
filename = "Box\neg1500_Box_p.xlsx";
neg1500_Box_p = readtable(filename);
neg1500_Box_p = table2array(neg1500_Box_p);

% Spike Computation
[U_x_pos, U_2_pos, ~, ~] = ProcessVel(pos1499_Spike_U,pos1500_Spike_U);
[U_x_neg, U_2_neg, ~, ~] = ProcessVel(neg1499_Spike_U,neg1500_Spike_U);

PosIntegral = U_x_pos(:,2) + U_2_pos + pos1500_Spike_p(:,4);
test = trapz(U_x_pos(:,1),pos1500_Spike_p(:,4));
PosIntegral = trapz(U_x_pos(:,1),PosIntegral);
NegIntegral = U_x_neg(:,2) + U_2_neg + neg1500_Spike_p(:,4);
NegIntegral = trapz(U_x_neg(:,1),NegIntegral);

Fd = -1 * (PosIntegral - NegIntegral);
Cd1 = Fd / (120/2);

% Spike Computation
[U_x_pos, U_2_pos, ~, ~] = ProcessVel(pos1499_Box_U,pos1500_Box_U);
[U_x_neg, U_2_neg, ~, ~] = ProcessVel(neg1499_Box_U,neg1500_Box_U);

PosIntegral = U_x_pos(:,2) + U_2_pos + pos1500_Box_p(:,4);
test = trapz(U_x_pos(:,1),pos1500_Box_p(:,4));
PosIntegral = trapz(U_x_pos(:,1),PosIntegral);
NegIntegral = U_x_neg(:,2) + U_2_neg + neg1500_Box_p(:,4);
NegIntegral = trapz(U_x_neg(:,1),NegIntegral);

Fd = -1 * (PosIntegral - NegIntegral);
Cd2 = Fd / (120/2);

% Functions
function [U_x, U_2,U1mag, U2mag]  = ProcessVel(U1,U2)
    U1mag = MagU(U1);
    U2mag = MagU(U2);
    U_x = zeros(size(U1mag,1),2);
    U_x(:,1) = U1(:,2);
    U_x(:,2) = U2mag - U1mag;
    U_2 = U2mag.^2;
end

function [mag] = MagU(U)
    mag = (U(:,4).^2 + U(:,5).^2).^(1/2);
end