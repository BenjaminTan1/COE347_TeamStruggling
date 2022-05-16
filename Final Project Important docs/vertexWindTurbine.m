%% COE 347 Final Project Block Mesh Dict
% Vertex Wind Turbine
clear all;
clc;

%% Set Dimensions
% Dimensions to Change
numSpikes = 0;                                % Number of spikes.
randSeed = 0;                                 % Randomness of spike distance.
yBox = 120;                                     % Length of box.
xBox = 800;                                     % Width of box.
xSpikes = 15;                                   % Dip of spikes.
xNeg = 2000;                                      % Front length of sim.
xPos = 2000;                                      % Back length of sim.
yTot = 1000;                                       % Upper edge of box.
% Meshing in Blocks
xMesh = 100;                                      % Mesh in x-direction of blocks.
yMesh = 100;                                     % Mesh in y-direction of blocks.
zMesh = 1;                                      % Mesh in z-direction of blocks.
expansionRatio = [2.0, 1.0, 1.0];               % Expansion ratio.

% Unchanging Dimensions
randStart = 1 - randSeed;
randStore = 1;
if numSpikes <= 0
    numSpikes = 0.5;
elseif numSpikes > 1
    randStore = randSeed*rand(numSpikes,1)+randStart;
end
ySpikes = yBox / (2*(numSpikes-1) + 2);
zDiff = 0.1;
meshing = [xMesh, yMesh, zMesh];

%% Generate Vertices
index = 1;
vertices = zeros(2*(16+(2*(numSpikes-1)+1)),3);
for z = 1:2
    % Inner Square and Spikes
    vertices(index,:) = [xBox/2, yBox/2, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [-xBox/2, yBox/2, (z-1)*zDiff];
    index = index + 1;
    inc = 1;
    for i = 1:2*(numSpikes-1)+1
        if mod(i,2)
            vertices(index,:) = [-xBox/2+(xSpikes*randStore(inc)), yBox/2-ySpikes*i, (z-1)*zDiff];
            inc = inc + 1;
        else
            vertices(index,:) = [-xBox/2, yBox/2-ySpikes*i, (z-1)*zDiff];
        end
        index = index + 1;
    end
    vertices(index,:) = [-xBox/2,-yBox/2, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [xBox/2,-yBox/2, (z-1)*zDiff];
    index = index + 1;

    % Outer edges of box
    % Top Edge
    vertices(index,:) = [xPos,yBox/2, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [xPos,yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [xBox/2,yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [-xBox/2,yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [-xNeg,yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [-xNeg,yBox/2, (z-1)*zDiff];
    index = index + 1;
    for i = 1:2*(numSpikes-1)+1
        vertices(index,:) = [-xNeg, yBox/2-ySpikes*i, (z-1)*zDiff];
        index = index + 1;
    end

    % Bot Edge
    vertices(index,:) = [-xNeg,-yBox/2, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [-xNeg,-yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [-xBox/2,-yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [xBox/2,-yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [xPos,-yTot, (z-1)*zDiff];
    index = index + 1;
    vertices(index,:) = [xPos,-yBox/2, (z-1)*zDiff];
    index = index + 1;
end

fileID=fopen('vertices.txt','w');
for i=1:size(vertices,1)
    fprintf(fileID,'(%e %e %e) \n',vertices(i,:));
end
fclose(fileID);

% Vertex Index of Spikes
spikes1 = [1, 3 + 2*(numSpikes-1)];
spikes2 = [10 + 2*(numSpikes-1), 11 + 2*(numSpikes-1)+1 + 2*(numSpikes-1)];

%% Generate Blocks
blocks = zeros(7+2*numSpikes,4);
blocks(1,:) = [spikes1(2) + 7, 1, spikes1(2) + 5, spikes1(2) + 6];
blocks(2,:) = [1, 0, spikes1(2) + 4, spikes1(2) + 5];
blocks(3,:) = [0, spikes1(2) + 2, spikes1(2) + 3, spikes1(2) + 4];
for i = 1:2*numSpikes
    blocks(3+i,:) = [spikes2(1)+(i-1)+1,spikes1(1)+(i-1)+1, spikes1(1)+(i-1), spikes2(1)+(i-1)];
end
blocks(4+2*numSpikes,:) = [spikes1(2)+1, spikes2(2)+5, spikes1(2)+2, 0];
blocks(5+2*numSpikes,:) = [spikes2(2)+1, spikes2(2)+2, spikes1(2), spikes2(2)];
blocks(6+2*numSpikes,:) = [spikes2(2)+2, spikes2(2)+3, spikes1(2)+1, spikes1(2)];
blocks(7+2*numSpikes,:) = [spikes2(2)+3, spikes2(2)+4, spikes2(2)+5, spikes1(2)+1];

blocks = [blocks, blocks+size(vertices,1)/2];

fileID=fopen('blocks.txt','w');
for i=1:size(blocks,1)
    fprintf(fileID,'hex (%d %d %d %d %d %d %d %d) (%d %d %d) simpleGrading ( %e %e %.1f ) \n',(blocks(i,:)),meshing, expansionRatio);
end
fclose(fileID);

%% Generate Faces
% Inlet
facesIn = zeros(2*numSpikes+2,4);
facesIn(1,:) = [spikes2(1)-1, spikes2(1)-1 + size(vertices,1)/2, spikes2(1) + size(vertices,1)/2, spikes2(1)];
for i = 1:2*numSpikes
    facesIn(1+i,:) = [spikes2(1)-1+i, spikes2(1)-1+i + size(vertices,1)/2, spikes2(1)+i + size(vertices,1)/2, spikes2(1)+i];
end
facesIn(end,:) = [spikes2(2), spikes2(2) + size(vertices,1)/2, spikes2(2)+1 + size(vertices,1)/2, spikes2(2)+1];

fileID=fopen('facesIn.txt','w');
fprintf(fileID, 'inlet\n{\n\ttype patch;\n\tfaces\n\t(\n');
for i=1:size(facesIn,1)
    fprintf(fileID,'\t\t(%d %d %d %d)\n', facesIn(i,:));
end
fprintf(fileID, '\t);\n}');
fclose(fileID);

% Outlet
facesOut = zeros(3,4);
facesOut(1,:) = [spikes1(2)+3, spikes1(2)+3+size(vertices,1)/2, spikes1(2)+2+size(vertices,1)/2, spikes1(2)+2];
facesOut(2,:) = [spikes1(2)+2, spikes1(2)+2+size(vertices,1)/2, spikes2(2)+5+size(vertices,1)/2, spikes2(2)+5];
facesOut(3,:) = [spikes2(2)+5, spikes2(2)+5+size(vertices,1)/2, spikes2(2)+4+size(vertices,1)/2, spikes2(2)+4];

fileID=fopen('facesOut.txt','w');
fprintf(fileID, 'outlet\n{\n\ttype patch;\n\tfaces\n\t(\n');
for i=1:size(facesOut,1)
    fprintf(fileID,'\t\t(%d %d %d %d)\n', facesOut(i,:));
end
fprintf(fileID, '\t);\n}');
fclose(fileID);

% Top
facesTop = zeros(3,4);
facesTop(1,:) = [spikes1(2)+6, spikes1(2)+6+size(vertices,1)/2, spikes1(2)+5+size(vertices,1)/2, spikes1(2)+5];
facesTop(2,:) = [spikes1(2)+5, spikes1(2)+5+size(vertices,1)/2, spikes1(2)+4+size(vertices,1)/2, spikes1(2)+4];
facesTop(3,:) = [spikes1(2)+4, spikes1(2)+4+size(vertices,1)/2, spikes1(2)+3+size(vertices,1)/2, spikes1(2)+3];

fileID=fopen('facesTop.txt','w');
fprintf(fileID, 'outlet\n{\n\ttype symmetryPlane;\n\tfaces\n\t(\n');
for i=1:size(facesTop,1)
    fprintf(fileID,'\t\t(%d %d %d %d)\n', facesTop(i,:));
end
fprintf(fileID, '\t);\n}');
fclose(fileID);

% Top
facesBot = zeros(3,4);
facesBot(1,:) = [spikes2(2)+1, spikes2(2)+1+size(vertices,1)/2, spikes2(2)+2+size(vertices,1)/2, spikes2(2)+2];
facesBot(2,:) = [spikes2(2)+2, spikes2(2)+2+size(vertices,1)/2, spikes2(2)+3+size(vertices,1)/2, spikes2(2)+3];
facesBot(3,:) = [spikes2(2)+3, spikes2(2)+3+size(vertices,1)/2, spikes2(2)+4+size(vertices,1)/2, spikes2(2)+4];

fileID=fopen('facesBot.txt','w');
fprintf(fileID, 'outlet\n{\n\ttype symmetryPlane;\n\tfaces\n\t(\n');
for i=1:size(facesBot,1)
    fprintf(fileID,'\t\t(%d %d %d %d)\n', facesBot(i,:));
end
fprintf(fileID, '\t);\n}');
fclose(fileID);

% Box
facesBox = zeros(2*numSpikes+3,4);
facesBox(1,:) = [0, size(vertices,1)/2, 1 + size(vertices,1)/2, 1];
for i = 1:2*numSpikes
    facesBox(1+i,:) = [spikes1(1)-1+i, spikes1(1)-1+i + size(vertices,1)/2, spikes1(1)+i + size(vertices,1)/2, spikes1(1)+i];
end
facesBox(end-1,:) = [spikes1(2), spikes1(2) + size(vertices,1)/2, spikes1(2)+1 + size(vertices,1)/2, spikes1(2)+1];
facesBox(end,:) = [spikes1(2)+1, spikes1(2)+1 + size(vertices,1)/2, size(vertices,1)/2, 0];

fileID=fopen('facesBox.txt','w');
fprintf(fileID, 'outlet\n{\n\ttype wall;\n\tfaces\n\t(\n');
for i=1:size(facesBox,1)
    fprintf(fileID,'\t\t(%d %d %d %d)\n', facesBox(i,:));
end
fprintf(fileID, '\t);\n}');
fclose(fileID);

%% Generate Figures
figure;
scatter(vertices(:,1),vertices(:,2), 10);
hold on;
for i = 1:(size(vertices,1)/2)
    if i > 2 && i <= 2 + 2*(numSpikes-1)+1 || i > 10 + 2*(numSpikes-1)+1 && i <= 10 + 2*(numSpikes-1)+1 + 2*(numSpikes-1)+1
    else
        text(vertices(i,1),vertices(i,2),num2str(i-1),'FontSize',8);
    end
end
% Blocks Z = 0
for i = 1:size(blocks,1)
    plot([vertices(blocks(i,1:4)+1,1);vertices(blocks(i,5)+1,1)], [vertices(blocks(i,1:4)+1,2);vertices(blocks(i,5)+1,2)]);
    hold on;
end
%hold off;

%figure;
scatter(vertices(:,1),vertices(:,2), 10);
hold on;
for i = 1:(size(vertices,1)/2)
    if i > 2 && i <= 2 + 2*(numSpikes-1)+1 || i > 10 + 2*(numSpikes-1)+1 && i <= 10 + 2*(numSpikes-1)+1 + 2*(numSpikes-1)+1
    else
        text(vertices(i,1),vertices(i,2),num2str(i+size(vertices,1)/2-1),'FontSize',8);
    end
end
% Blocks Z = z0
for i = 1:size(blocks,1)
    plot([vertices(blocks(i,5:8)+1,1);vertices(blocks(i,5)+1,1)], [vertices(blocks(i,5:8)+1,2);vertices(blocks(i,5)+1,2)]);
    hold on;
end

