% Vertex Wind Turbine
clc;

% Dimensions to Change
numSpikes = 3;
if numSpikes <= 0
    numSpikes = 0.5;
end
yBox = 1;                                   % Length of box.
xBox = 1;                                   % Width of box.
xSpikes = 0.1;                             % Dip of spikes.
r0 = 1.5;                                   % Boundary layer radius.
xNeg = 5;                                   % Front length of box.
xPos = 12;                                  % Back length of box.
yTot = 7;                                   % Upper edge of box.

% Set Dimensions
ySpikes = yBox / (2*(numSpikes-1) + 2);
zDiff = 0.1;
r = r0 * max(yBox,xBox);

vertices = zeros(2*(24+(2*(numSpikes-1)+1)),3);
for z = 1:2
    % Inner Square and Spikes
    vertices(1+(z-1)*size(vertices,1)/2,:) = [xBox/2, yBox/2, (z-1)*zDiff];
    vertices(2+(z-1)*size(vertices,1)/2,:) = [-xBox/2, yBox/2, (z-1)*zDiff];
    for i = 1:2*(numSpikes-1)+1
        if mod(i,2)
            vertices((i+2)+(z-1)*size(vertices,1)/2,:) = [-xBox/2+xSpikes, yBox/2-ySpikes*i, (z-1)*zDiff];
        else
            vertices((i+2)+(z-1)*size(vertices,1)/2,:) = [-xBox/2, yBox/2-ySpikes*i, (z-1)*zDiff];
        end
    end
    vertices(3+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xBox/2,-yBox/2, (z-1)*zDiff];
    vertices(4+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xBox/2,-yBox/2, (z-1)*zDiff];

    % Radial Circle
    thetaX = acos(xBox/r);
    thetaY = acos(yBox/r);
    
    for quad = 1:4

        if quad <= 2
            signY = 1;
        else
            signY = -1;
        end

        if quad >= 2 && quad <= 3
            signX = -1;
        else
            signX = 1;
        end

        boolOdd = mod(quad,2);
        boolEven = ~mod(quad,2);
        vertices(boolOdd+5+2*(quad-1)+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [signX*r*sin(thetaY), signY * yBox/2, (z-1)*zDiff];
        vertices(boolEven+5+2*(quad-1)+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [signX*xBox/2, signY * r*sin(thetaX), (z-1)*zDiff];
    end

    % Outer edges of box
    % Top Edge
    vertices(13+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xPos,yBox/2, (z-1)*zDiff];
    vertices(14+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xPos,yTot, (z-1)*zDiff];
    vertices(15+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xBox/2,yTot, (z-1)*zDiff];
    vertices(16+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xBox/2,yTot, (z-1)*zDiff];
    vertices(17+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xNeg,yTot, (z-1)*zDiff];
    vertices(18+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xNeg,yBox/2, (z-1)*zDiff];
    % Bot Edge
    vertices(19+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xNeg,-yBox/2, (z-1)*zDiff];
    vertices(20+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xNeg,-yTot, (z-1)*zDiff];
    vertices(21+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [-xBox/2,-yTot, (z-1)*zDiff];
    vertices(22+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xBox/2,-yTot, (z-1)*zDiff];
    vertices(23+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xPos,-yTot, (z-1)*zDiff];
    vertices(24+2*(numSpikes-1)+1+(z-1)*size(vertices,1)/2,:) = [xPos,-yBox/2, (z-1)*zDiff];
end

figure;
scatter(vertices(:,1),vertices(:,2));
