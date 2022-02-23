%% OpenFOAM Wall Shear Stress Calculations
clear;

%% Initializing File Names
% Replace with your own file path for this data.
% Square Data
file_name1 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\square_data\Data1_RE10";
file_name2 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\square_data\Data2_RE100";
file_name3 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\square_data\Data3_RE250";
file_name4 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\square_data\Data4_RE500";

% Shallow Data
file_name5 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\shallow_data\shallowRe10_data";
file_name6 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\shallow_data\shallowRe50_data";
file_name7 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\shallow_data\shallowRe100_data";
file_name8 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\shallow_data\shallowRe250_data";
file_name9 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\shallow_data\shallowRe500_data";

% Tall Data
file_name10 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\tall_data\tallRe10_data";
file_name11 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\tall_data\tallRe50_data";
file_name12 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\tall_data\tallRe100_data";
file_name13 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\tall_data\tallRe250_data";
file_name14 = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\tall_data\tallRe500_data";

%% Function Calls
% Square
[u_x1, u_y1, ~, dp_y1, F1, ~] = GridXY(file_name1, 81, 0.1, 0.1);
[u_x2, u_y2, ~, dp_y2, F2, ~] = GridXY(file_name2, 81, 0.1, 0.1);
[u_x3, u_y3, ~, dp_y3, F3, ~] = GridXY(file_name3, 81, 0.1, 0.1);
[u_x4, u_y4, ~, dp_y4, F4, ~] = GridXY(file_name4, 100, 0.1, 0.1);
% [x_x, p_x, dp_x, dUx] = DifferentiateGridX(u_x);
% [x_y, p_y, dp_y, dUy] = DifferentiateGridY(u_x);

% Shallow
[u_x5, u_y5, ~, dp_y5, F5, ~] = GridXY(file_name5, 81, 0.2, 0.1);
[u_x6, u_y6, ~, dp_y6, F6, ~] = GridXY(file_name6, 81, 0.2, 0.1);
[u_x7, u_y7, ~, dp_y7, F7, ~] = GridXY(file_name7, 81, 0.2, 0.1);
[u_x8, u_y8, ~, dp_y8, F8, ~] = GridXY(file_name8, 81, 0.2, 0.1);
[u_x9, u_y9, ~, dp_y9, F9, ~] = GridXY(file_name9, 81, 0.2, 0.1);

% Tall
[u_x10, u_y10, ~, dp_y10, F10, ~] = GridXY(file_name10, 81, 0.1, 0.2);
[u_x11, u_y11, ~, dp_y11, F11, ~] = GridXY(file_name11, 81, 0.1, 0.2);
[u_x12, u_y12, ~, dp_y12, F12, ~] = GridXY(file_name12, 81, 0.1, 0.2);
[u_x13, u_y13, ~, dp_y13, F13, ~] = GridXY(file_name13, 81, 0.1, 0.2);
[u_x14, u_y14, ~, dp_y14, F14, ~] = GridXY(file_name14, 81, 0.1, 0.2);

%% Generates figures.
% Define x and y values.
x81 = linspace(0, 0.1, 81);
y81 = linspace(0, 0.1, 81);
x81s = linspace(0, 0.2, 81);
y81s = linspace(0, 0.1, 81);
x100 = linspace(0, 0.1, 100);
y100 = linspace(0, 0.1, 100);
x81t = linspace(0, 0.1, 81);
y81t = linspace(0, 0.2, 81);

%% Square
% % Contour Plots of Velocity
% figure;
% contour(x100, y100, u_x4);
% title('Plot of u');
% xlabel('x-distance (m)');
% ylabel('y-distace (m)');
% 
% figure;
% contour(x100, y100, u_y4);
% title('Plot of v');
% xlabel('x-distance (m)');
% ylabel('y-distance (m)');

% Shear Force Graphs
figure;
plot(x81, dp_y1);
hold on;
plot(x81, dp_y2);
hold on;
plot(x81, dp_y3);
hold on;
plot(x100, dp_y4);
hold on;
title('Plot of Shear Stress on Boundary');
xlabel('x-distance (m)');
ylabel('Tal (a.u.)');
hold off;
legend('Re = 10','Re = 100', 'Re = 250', 'Re = 500');

% Non-Dimensional Force Graphs
F = [F1, F2, F3, F4];
Re = [10, 100, 250, 500];

figure;
for i = 1:size(Re, 2)
    scatter(Re(i), F(i));
    hold on;
end
plot(Re, F);
title('Plot of Non-Dimensional Force Over Reynolds Number');
legend('Re = 10', 'Re = 100', 'Re = 250', 'Re = 500');
xlabel('Reynolds Number');
ylabel('Intensity (a.u.)');

%% Shallow
F = [F5, F6, F7, F8, F9];
Re = [10, 50, 100, 250, 500];

% Shear Stress Graphs
figure;
plot(x81s, dp_y5);
hold on;
plot(x81s, dp_y6);
hold on;
plot(x81s, dp_y7);
hold on;
plot(x81s, dp_y8);
hold on;
plot(x81s, dp_y9);
hold on;
title('Plot of Shear Stress on Boundary (Shallow 0.2 x 0.1 m)');
xlabel('x-distance (m)');
ylabel('Tal (a.u.)');
hold off;
legend('Re = 10', 'Re = 50', 'Re = 100', 'Re = 250', 'Re = 500');

% Non-Dimensional Force Graphs
figure;
for i = 1:size(Re, 2)
    scatter(Re(i), F(i));
    hold on;
end
plot(Re, F);
title('Plot of Non-Dimensional Force Over Reynolds Number (Shallow 0.2 x 0.1 m)');
legend('Re = 10', 'Re = 50', 'Re = 100', 'Re = 250', 'Re = 500');
xlabel('Reynolds Number');
ylabel('Intensity (a.u.)');

%% Tall
F = [F10, F11, F12, F13, F14];
Re = [10, 50, 100, 250, 500];

figure;
plot(x81, dp_y10);
hold on;
plot(x81, dp_y11);
hold on;
plot(x81, dp_y12);
hold on;
plot(x81, dp_y13);
hold on;
plot(x81, dp_y14);
hold on;
title('Plot of Shear Stress on Boundary (Tall 0.1 x 0.2 m)');
xlabel('x-distance (m)');
ylabel('Tal (a.u.)');
hold off;
legend('Re = 10', 'Re = 50', 'Re = 100', 'Re = 250', 'Re = 500');

% Non-Dimensional Force Graphs
figure;
for i = 1:size(Re, 2)
    scatter(Re(i), F(i));
    hold on;
end
plot(Re, F);
title('Plot of Non-Dimensional Force Over Reynolds Number (Tall 0.1 x 0.2 m)');
legend('Re = 10', 'Re = 50', 'Re = 100', 'Re = 250', 'Re = 500');
xlabel('Reynolds Number');
ylabel('Intensity (a.u.)');

% Contour Plots of Velocity
figure;
contour(x81, y81, u_y1);
title('Plot of uy for Square Re = 10');
xlabel('x-distance (m)');
ylabel('y-distace (m)');

figure;
contour(x81s, y81s, u_y5);
title('Plot of uy for Shallow Re = 10');
xlabel('x-distance (m)');
ylabel('y-distace (m)');

figure;
contour(x81t, y81t, u_y10);
title('Plot of uy for Tall Re = 10');
xlabel('x-distance (m)');
ylabel('y-distace (m)');

% % Plotting the polynomial generated by the DifferentiateGridX function.
% for i = 1:21
%     figure;
%     plot(x_x, polyval(p_x(:,i), x_x));
%     strTitle = ["Polynomial for", 0.05 * (i-1)];
%     title(strTitle);
% end
% 
% % Plotting the polynomial generated by the DifferentiateGridX function.
% for i = 1:21
%     figure;
%     plot(x_y, polyval(p_y(:,i), x_y));
%     strTitle = ["Polynomial for", 0.05 * (i-1), " m"];
%     title(strTitle);
% end

%% Functions
% Tabulates the values of U_x and U_y into a matrix containing the
% component of velocity for each mesh row in a corresponding row.
% Calls DifferentiateGridY for returning the polynomial fit and derivative
% of the polynomial fit.
function [u_xf, u_yf, p_x_yf, duf_y, F, Tf] = GridXY(file_name, meshsize, L, H)

    T = readtable(file_name);
    T = table2array(T);
    
    u_x = zeros(meshsize, meshsize);
    u_y = zeros(meshsize, meshsize);
    du_y = zeros(meshsize, meshsize);
    
    for j = 1:meshsize
        for i = 1:meshsize
            u_x(j, i) = T((meshsize*(j-1))+i,1);
            u_y(j, i) = T((meshsize*(j-1))+i,2);
        end
    end
    
    % Generate polynomial for u_x going down the y-direction.
    [y, p_x_yf, dp_x_y, ~] = DifferentiateGridY(u_x, H);
    
    % Evaluate the derivative of the polynomial for u_x going down the
    % y-direction.
    for i = 1:meshsize
        du_y(:, i) = polyval(dp_x_y(:,i), H);
    end
    
    % Integrate the mesh of the derivative using trapezoidal method.
    x_mesh = linspace(0, L, meshsize);
    F = trapz(x_mesh, du_y(meshsize,:));
    
    % Returns output.
    u_xf = u_x;
    u_yf = u_y;
    duf_y = du_y(meshsize,:);
    Tf = T;
end

% Creates a polynomial fit of a matrix of values going through its rows.
function [xf, pf, dpf, dUf] = DifferentiateGridX(U)

    % Initialize variables.
    [m, n] = size(U);
    p = zeros(m, n);
    dU = zeros(m, n);
    dp = zeros(m-1, n);
    x = linspace(0, 0.1, n);
    
    % Create polynomial equations.
    for i = 1:m
        p(:,i) = polyfit(x, U(i,:), m-1);
        dp(:,i) = polyder(p(:,i));
    end
    
    % Returns output
    xf = x;
    pf = p;
    dpf = dp;
    dUf = 1;
end

% Creates a polynomial fit of a matrix of values going down its columns.
% Finds 
function [yf, pf, dpf, dUf] = DifferentiateGridY(U, H)

    % Initialize variables.
    [m, n] = size(U);
    p = zeros(6, n);
    dU = zeros(m, n);
    dp = zeros(5, n);
    y = linspace(0, H, m);
    
    % Create polynomial equations.
    for j = 1:n
        p(:,j) = polyfit(y, U(:,j), 5);
        dp(:,j) = polyder(p(:,j));
    end
    
    % Returns output.
    yf = y;
    pf = p;
    dpf = dp;
    dUf = 1;
end