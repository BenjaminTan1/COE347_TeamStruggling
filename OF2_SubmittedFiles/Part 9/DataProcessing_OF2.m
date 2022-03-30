%% COE 347 Data Processing OF2
clear all;

%% Extra File Data.
% Refined
file_name1 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\line_U_div4";
file_name2 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\line_U_div2";
file_name3 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\line_U_3div4";

T1 = readtable(file_name1);
T2 = readtable(file_name2);
T2(2:2:end,:) = [];
T3 = readtable(file_name3);

T1 = table2array(T1);
T2 = table2array(T2);
T3 = table2array(T3);

% Medium
file_name4 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\Medium\line1_U";
file_name5 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\Medium\line2_U";
file_name6 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\Medium\line3_U";

T4 = readtable(file_name4);
T5 = readtable(file_name5);
T5(3:2:end,:) = [];
T6 = readtable(file_name6);

T4 = table2array(T4);
T5 = table2array(T5);
T6 = table2array(T6);

% Coarse
file_name7 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\Coarse\line1_U";
file_name8 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\Coarse\line2_U";
file_name9 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\OF2\Coarse\line3_U";

T7 = readtable(file_name7);
T8 = readtable(file_name8);
T8(3:2:end,:) = [];
T9 = readtable(file_name9);

T7 = table2array(T7);
T8 = table2array(T8);
T9 = table2array(T9);

%% Data Processing.
% Refined
[~, ~, ~, rPI4, u_xPI4, u_yPI4, u_zPI4, u_rPI4, u_thetaPI4, thetaPI4] = DataProcessing(T1);
[e_rrPI4, e_rthetaPI4] = Strain(rPI4, u_rPI4, u_thetaPI4);

[~, ~, ~, rPI2, u_xPI2, u_yPI2, u_zPI2, u_rPI2, u_thetaPI2, thetaPI2] = DataProcessing(T2);
[e_rrPI2, e_rthetaPI2] = Strain(rPI2, u_rPI2, u_thetaPI2);

[~, ~, ~, r3PI4, u_x3PI4, u_y3PI4, u_z3PI4, u_r3PI4, u_theta3PI4, theta3PI4] = DataProcessing(T3);
[e_rr3PI4, e_rtheta3PI4] = Strain(r3PI4, u_r3PI4, u_theta3PI4);

u_thetaArr = [u_thetaPI4 u_thetaPI2 u_theta3PI4];
r_Arr = [rPI4 rPI2 r3PI4];
thetaMesh = [thetaPI4, thetaPI2, theta3PI4];
[C, F, du_theta] = DragCoefficient(r_Arr, u_thetaArr, thetaMesh, 1, 2, pi);

% Medium
[~, ~, ~, rPI4m, u_xPI4m, u_yPI4m, u_zPI4m, u_rPI4m, u_thetaPI4m, thetaPI4m] = DataProcessing(T4);
[e_rrPI4m, e_rthetaPI4m] = Strain(rPI4m, u_rPI4m, u_thetaPI4m);

[~, ~, ~, rPI2m, u_xPI2m, u_yPI2m, u_zPI2m, u_rPI2m, u_thetaPI2m, thetaPI2m] = DataProcessing(T5);
[e_rrPI2m, e_rthetaPI2m] = Strain(rPI2m, u_rPI2m, u_thetaPI2m);

[~, ~, ~, r3PI4m, u_x3PI4m, u_y3PI4m, u_z3PI4m, u_r3PI4m, u_theta3PI4m, theta3PI4m] = DataProcessing(T6);
[e_rr3PI4m, e_rtheta3PI4m] = Strain(r3PI4m, u_r3PI4m, u_theta3PI4m);

u_thetaArrm = [u_thetaPI4m u_thetaPI2m u_theta3PI4m];
r_Arr = [rPI4m rPI2m r3PI4m];
thetaMesh = [thetaPI4m, thetaPI2m, theta3PI4m];
[Cm, Fm, du_thetam] = DragCoefficient(r_Arr, u_thetaArrm, thetaMesh, 1, 2, pi);

% Coarse
[~, ~, ~, rPI4c, u_xPI4c, u_yPI4c, u_zPI4c, u_rPI4c, u_thetaPI4c, thetaPI4c] = DataProcessing(T7);
[e_rrPI4c, e_rthetaPI4c] = Strain(rPI4c, u_rPI4c, u_thetaPI4c);

[~, ~, ~, rPI2c, u_xPI2c, u_yPI2c, u_zPI2c, u_rPI2c, u_thetaPI2c, thetaPI2c] = DataProcessing(T8);
[e_rrPI2c, e_rthetaPI2c] = Strain(rPI2c, u_rPI2c, u_thetaPI2c);

[~, ~, ~, r3PI4c, u_x3PI4c, u_y3PI4c, u_z3PI4c, u_r3PI4c, u_theta3PI4c, theta3PI4c] = DataProcessing(T9);
[e_rr3PI4c, e_rtheta3PI4c] = Strain(r3PI4c, u_r3PI4c, u_theta3PI4c);

u_thetaArrc = [u_thetaPI4c u_thetaPI2c u_theta3PI4c];
r_Arr = [rPI4c rPI2c r3PI4c];
thetaMesh = [thetaPI4c, thetaPI2c, theta3PI4c];
[Cc, Fc, du_thetac] = DragCoefficient(r_Arr, u_thetaArrc, thetaMesh, 1, 2, pi);

%% Plot Figures.
% Velocity r
figure;
plot(rPI4, u_rPI4);
hold on;
plot(rPI2, u_rPI2);
hold on;
plot(-r3PI4, u_r3PI4);
title('Refined Velocity r')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4m, u_rPI4m);
hold on;
plot(rPI2m, u_rPI2m);
hold on;
plot(-r3PI4m, u_r3PI4m);
title('Medium Velocity r')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4c, u_rPI4c);
hold on;
plot(rPI2c, u_rPI2c);
hold on;
plot(-r3PI4c, u_r3PI4c);
title('Coarse Velocity r')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

% Plot Near Wall
figure;
plot(rPI4, u_rPI4);
hold on;
plot(rPI2, u_rPI2);
hold on;
plot(-r3PI4, u_r3PI4);
hold on;
plot(rPI4m, u_rPI4m);
hold on;
plot(rPI2m, u_rPI2m);
hold on;
plot(-r3PI4m, u_r3PI4m);
hold on;
plot(rPI4c, u_rPI4c);
hold on;
plot(rPI2c, u_rPI2c);
hold on;
plot(-r3PI4c, u_r3PI4c);
title('Velocity r')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('refined Pi/4','refined Pi/2','refined 3Pi/4', 'medium Pi/4', 'medium Pi/2', 'medium 3Pi/4', 'coarse Pi/4', 'coarse Pi/2', 'coarse 3Pi/4');
xlim([0.5 0.6]);

% Velocity theta
figure;
plot(rPI4, u_thetaPI4);
hold on;
plot(rPI2, u_thetaPI2);
hold on;
plot(-r3PI4, u_theta3PI4);
title('Refined Velocity theta')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4m, u_thetaPI4m);
hold on;
plot(rPI2m, u_thetaPI2m);
hold on;
plot(-r3PI4m, u_theta3PI4m);
title('Medium Velocity theta')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4c, u_thetaPI4c);
hold on;
plot(rPI2c, u_thetaPI2c);
hold on;
plot(-r3PI4c, u_theta3PI4c);
title('Coarse Velocity theta')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

% Plot Near Wall
figure;
plot(rPI4, u_thetaPI4);
hold on;
plot(rPI2, u_thetaPI2);
hold on;
plot(-r3PI4, u_theta3PI4);
hold on;
plot(rPI4m, u_thetaPI4m);
hold on;
plot(rPI2m, u_thetaPI2m);
hold on;
plot(-r3PI4m, u_theta3PI4m);
hold on;
plot(rPI4c, u_thetaPI4c);
hold on;
plot(rPI2c, u_thetaPI2c);
hold on;
plot(-r3PI4c, u_theta3PI4c);
title('Velocity theta')
xlabel('Position (m)')
ylabel('Velocity (m/s)');
legend('refined Pi/4','refined Pi/2','refined 3Pi/4', 'medium Pi/4', 'medium Pi/2', 'medium 3Pi/4', 'coarse Pi/4', 'coarse Pi/2', 'coarse 3Pi/4');
xlim([0.5 0.6]);

% Strain r
figure;
plot(rPI4, e_rrPI4);
hold on;
plot(rPI2, e_rrPI2);
hold on;
plot(-r3PI4, e_rr3PI4);
title('Refined Strain Tensor r,r')
xlabel('Position (m)')
ylabel('Strain');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4m, e_rrPI4m);
hold on;
plot(rPI2m, e_rrPI2m);
hold on;
plot(-r3PI4m, e_rr3PI4m);
title('Medium Strain Tensor r,r')
xlabel('Position (m)')
ylabel('Strain');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4c, e_rrPI4c);
hold on;
plot(rPI2c, e_rrPI2c);
hold on;
plot(-r3PI4c, e_rr3PI4c);
title('Coarse Strain Tensor r,r')
xlabel('Position (m)')
ylabel('Strain');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4, e_rrPI4);
hold on;
plot(rPI2, e_rrPI2);
hold on;
plot(-r3PI4, e_rr3PI4);
hold on;
plot(rPI4m, e_rrPI4m);
hold on;
plot(rPI2m, e_rrPI2m);
hold on;
plot(-r3PI4m, e_rr3PI4m);
hold on;
plot(rPI4c, e_rrPI4c);
hold on;
plot(rPI2c, e_rrPI2c);
hold on;
plot(-r3PI4c, e_rr3PI4c);
title('Strain Tensor r,r')
xlabel('Position (m)')
ylabel('Strain');
legend('refined Pi/4','refined Pi/2','refined 3Pi/4', 'coarse Pi/4','coarse Pi/2','coarse 3Pi/4', 'medium Pi/4','medium Pi/2','medium 3Pi/4');
xlim([0.5 1]);

% Strain theta
figure;
plot(rPI4, e_rthetaPI4);
hold on;
plot(rPI2, e_rthetaPI2);
hold on;
plot(-r3PI4, e_rtheta3PI4);
title('Refined Strain Tensor r,theta')
xlabel('Position (m)')
ylabel('Strain');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4c, e_rthetaPI4c);
hold on;
plot(rPI2c, e_rthetaPI2c);
hold on;
plot(-r3PI4c, e_rtheta3PI4c);
title('Coarse Strain Tensor r,theta')
xlabel('Position (m)')
ylabel('Strain');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4m, e_rthetaPI4m);
hold on;
plot(rPI2m, e_rthetaPI2m);
hold on;
plot(-r3PI4m, e_rtheta3PI4m);
title('Medium Strain Tensor r,theta')
xlabel('Position (m)')
ylabel('Strain');
legend('Pi/4','Pi/2','3Pi/4');
xlim([0.5 1]);

figure;
plot(rPI4, e_rthetaPI4);
hold on;
plot(rPI2, e_rthetaPI2);
hold on;
plot(-r3PI4, e_rtheta3PI4);
hold on;
plot(rPI4c, e_rthetaPI4c);
hold on;
plot(rPI2c, e_rthetaPI2c);
hold on;
plot(-r3PI4c, e_rtheta3PI4c);
hold on;
plot(rPI4m, e_rthetaPI4m);
hold on;
plot(rPI2m, e_rthetaPI2m);
hold on;
plot(-r3PI4m, e_rtheta3PI4m);
title('Strain Tensor r,theta')
xlabel('Position (m)')
ylabel('Strain');
legend('refined Pi/4','refined Pi/2','refined 3Pi/4', 'coarse Pi/4','coarse Pi/2','coarse 3Pi/4', 'medium Pi/4','medium Pi/2','medium 3Pi/4');
xlim([0.5 1]);

% Strain Sigma
e_rrPI4Compiled = [e_rrPI4m(2), e_rrPI4(2), e_rrPI4c(2)];
e_rrPI2Compiled = [e_rrPI2m(2), e_rrPI2(2), e_rrPI2c(2)];
e_rr3PI4Compiled = [e_rr3PI4m(2), e_rr3PI4(2), e_rr3PI4c(2)];
Sigma = [ComputeSigma(rPI4m), ComputeSigma(rPI4), ComputeSigma(rPI4c)];

figure;
plot(Sigma, e_rrPI4Compiled);
hold on;
plot(Sigma, e_rrPI2Compiled);
hold on;
plot(Sigma, e_rr3PI4Compiled);
title('Strain rr vs. Mesh Size');
xlabel('Sigma')
ylabel('Strain Tensor rr');
legend('Pi/4','Pi/2','3Pi/4');

e_rthetaPI4Compiled = [e_rthetaPI4m(2), e_rthetaPI4(2), e_rthetaPI4c(2)];
e_rthetaPI2Compiled = [e_rthetaPI2m(2), e_rthetaPI2(2), e_rthetaPI2c(2)];
e_rtheta3PI4Compiled = [e_rtheta3PI4m(2), e_rtheta3PI4(2), e_rtheta3PI4c(2)];
Sigma = [ComputeSigma(rPI4m), ComputeSigma(rPI4), ComputeSigma(rPI4c)];

figure;
plot(Sigma, e_rthetaPI4Compiled);
hold on;
plot(Sigma, e_rthetaPI2Compiled);
hold on;
plot(Sigma, e_rtheta3PI4Compiled);
title('Strain r,theta vs. Mesh Size');
xlabel('Sigma')
ylabel('Strain Tensor rr');
legend('Pi/4','Pi/2','3Pi/4');

file_nameL1 = 'L_Refined';
file_nameL2 = 'L_Medium';
file_nameL3 = 'L_Coarse';

% Recirculation Region

[rr, Lr] = ComputeL(file_nameL1);
[rm, Lm] = ComputeL(file_nameL2);
[rc, Lc] = ComputeL(file_nameL3);

LrF = 1.6632/pi;
LmF = 1.6209/pi;
LcF = 1.6203/pi;

L = [LrF LmF LcF];
mesh = [35, 30, 25];
figure;
scatter(mesh, L);
title('L/D Over Mesh Size');
xlabel('Mesh Size');
ylabel('L/D');

% Force

figure;
angles = [pi/4, pi/2, 3*pi/4];
du_graph = [du_theta];
du_graphm = [du_thetam];
du_graphc = [du_thetac];
plot(angles, du_graph)
hold on;
plot(angles, du_graphm)
hold on;
plot(angles, du_graphc)
title('Plot of du, theta, r at Different Mesh Sizes over Angle');
xlabel('Angle (rad)');
ylabel('du, theta, r');
legend('refined','medium','coarse');

%% Functions
% Returns components of velocity in both cartesian and radial coordinates.
function [x, y, z, r, u_x, u_y, u_z, u_rF, u_thetaF, thetaF] = DataProcessing(T)

    [m, ~] = size(T);
    
    % Return coordinates of position and coordinates of velocity.
    x = T(:,1);
    y = T(:,2);
    z = T(:,3);

    u_x = T(:,4);
    u_y = T(:,5);
    u_z = T(:,6);

    % Conversion into polar coordinates.
    r = zeros(m,1);
    u_r = zeros(m,1);
    u_theta = zeros(m,1);

    for i = 1:m
        theta = atan(y(i)/x(i));
        
        if theta == pi/2
            r(i) = y(i) / sin(theta);
        else
            r(i) = x(i) / cos(theta);
        end

        u_r(i) = cos(theta) * u_x(i) + sin(theta) * u_y(i);
        u_theta(i) = -sin(theta) * u_x(i) + cos(theta) * u_y(i);
    end

    if theta < 0
        theta = theta + pi;
    end

    % Return statements.
    u_rF = u_r;
    u_thetaF = u_theta;
    thetaF = theta;

end

function [e_rr, e_rtheta] = Strain(r, u_r, u_theta)

    % Evaluation for e_rr.
    p_rr = polyfit(r, u_r, 12);
    dp_rr = polyder(p_rr);

    % Evaluation for e_rtheta.
    p_rtheta = polyfit(r, u_theta./r, 12);
    dp_rtheta = polyder(p_rtheta);
    temp = polyval(dp_rtheta, r);
    temp = (r/2).*temp;

    % Return statements.
    e_rr = polyval(dp_rr,r);
    e_rtheta = temp;

end

function [C, F, du_theta] = DragCoefficient(r, u_thetaArr, thetaMesh, rho, U, D)
    
    [~, n] = size(u_thetaArr);
    
    % Set degree of poylnomial for interpolation.
    H = 12;
    p = zeros(H+1,n);
    dp = zeros(H,n);
    du = zeros(1,n);

    for i = 1:n 
        p(:,i) = polyfit(r(:,i), u_thetaArr(:,i), H);
        dp(:,i) = polyder(p(:,i));
        du(1, i) = polyval(dp(:,i), r(1,i));
    end
    
    F = 2*trapz(thetaMesh, du);
    du_theta = du;
    C = 2 * F / (rho * U^2 * D);

end

function [Sigma] = ComputeSigma(r)
    Sigma = abs((r(2) - r(1))*(r(3)-r(2)));
end

function [r, L] = ComputeL(file_name)

    % Degree of polynomial interpolation
    H = 4;

    T = readtable(file_name);
    T = table2array(T);

    p = polyfit(T(:,1), T(:,2),H);
    r = roots(p);

    L = r;
end