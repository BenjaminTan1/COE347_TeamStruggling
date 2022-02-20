%% OpenFOAM Wall Shear Stress Calculations
clear;

% Initializing File Names
file_name = "C:\Users\Benjamin Tan\Desktop\Courses\COE347\OpenFOAM\Data1_RE10";

% Function Calls
[u_x, u_y, T] = GridXY(file_name, 21);
[x, p, dUf] = DifferentiateGridY(u_x);

% Contour Plots of Velocity
figure;
contour(u_x);
figure;
contour(u_y);

for i = 1:21
    figure;
    plot(x, polyval(p(:,i), x));
    strTitle = ["Polynomial for", 0.05 * (i-1)];
    title(strTitle);
end

%% Functions
% Tabulates the values of U_x and U_y into a matrix containing the
% component of velocity for each mesh row in a corresponding row.
function [u_xf, u_yf, Tf] = GridXY(file_name, meshsize)

    T = readtable(file_name);
    T = table2array(T);
    
    u_x = zeros(meshsize, meshsize);
    u_y = zeros(meshsize, meshsize);
    
    for j = 1:meshsize
        for i = 1:meshsize
            u_x(j, i) = T((meshsize*(j-1))+i,1);
            u_y(j, i) = T((meshsize*(j-1))+i,2);
        end
    end
    
    u_xf = u_x;
    u_yf = u_y;
    
    Tf = T;
    
end

function [xf, pf, dUf] = DifferentiateGridY(U)

    [m, n] = size(U);
    p = zeros(m, n);
    x = linspace(0, 0.1, m);
    
    for j = 1:n
        p(:,j) = polyfit(x, U(:,j), n-1);
    end
    
    pf = p;
    dUf = 1;
    xf = x;
end

function [xf, pf, dUf] = DifferentiateGridX(U)

    [m, n] = size(U);
    p = zeros(m, n);
    x = linspace(0, 0.1, n);
    
    for i = 1:m
        p(:,i) = polyfit(x, U(i,:), m-1);
    end
    
    pf = p;
    dUf = 1;
    xf = x;
end