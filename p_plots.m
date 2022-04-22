clear
close all
clc

fig = 0;

%% Figure 1: Top Wall
meshAtopWall = importdata('meshA_topWall_p_T.xy');
meshBtopWall = importdata('meshB_topWall_p_T.xy');
meshCtopWall = importdata('meshC_topWall_p_T.xy');
meshDtopWall = importdata('meshD_topWall_p_T.xy');

fig = fig+1;
figure(fig)
hold on

plot(meshAtopWall(:,1),meshAtopWall(:,4),'blue');
plot(meshBtopWall(:,1),meshBtopWall(:,4),'red');
plot(meshCtopWall(:,1),meshCtopWall(:,4),'green');
plot(meshDtopWall(:,1),meshDtopWall(:,4),'magenta');
xlabel('x');
ylabel('Pressure, {\it p}');
legend ('Mesh A','Mesh B', 'Mesh C', 'Mesh D');
title('Spatial Distribution of Pressure at t = 4 on the Top Wall');
hold off

%% Figure 2: Step Horizontal
meshAstepHorizontal = importdata('meshA_stepHorizontal_p_T.xy');
meshBstepHorizontal = importdata('meshB_stepHorizontal_p_T.xy');
meshCstepHorizontal = importdata('meshC_stepHorizontal_p_T.xy');
meshDstepHorizontal = importdata('meshD_stepHorizontal_p_T.xy');

fig = fig+1;
figure(fig)
hold on

plot(meshAstepHorizontal(:,1),meshAstepHorizontal(:,4),'blue');
plot(meshBstepHorizontal(:,1),meshBstepHorizontal(:,4),'red');
plot(meshCstepHorizontal(:,1),meshCstepHorizontal(:,4),'green');
plot(meshDstepHorizontal(:,1),meshDstepHorizontal(:,4),'magenta');
xlabel('x');
ylabel('Pressure, {\it p}');
legend ('Mesh A','Mesh B', 'Mesh C', 'Mesh D');
title('Spatial Distribution of Pressure at t = 4 on the Step Horizontal');
hold off

%% Figure 3: Step Face
meshAstepFace = importdata('meshA_stepFace_p_T.xy');
meshBstepFace = importdata('meshB_stepFace_p_T.xy');
meshCstepFace = importdata('meshC_stepFace_p_T.xy');
meshDstepFace = importdata('meshD_stepFace_p_T.xy');

fig = fig+1;
figure(fig)
hold on

plot(meshAstepFace(:,2),meshAstepFace(:,4),'blue');
plot(meshBstepFace(:,2),meshBstepFace(:,4),'red');
plot(meshCstepFace(:,2),meshCstepFace(:,4),'green');
plot(meshDstepFace(:,2),meshDstepFace(:,4),'magenta');
xlabel('y');
ylabel('Pressure, {\it p}');
legend ('Mesh A','Mesh B', 'Mesh C', 'Mesh D');
title('Spatial Distribution of Pressure at t = 4 on the Step Face');
hold off

%% Figure 4: Entry Wall
meshAentryWall = importdata('meshA_entryWall_p_T.xy');
meshBentryWall = importdata('meshB_entryWall_p_T.xy');
meshCentryWall = importdata('meshC_entryWall_p_T.xy');
meshDentryWall = importdata('meshD_entryWall_p_T.xy');

fig = fig+1;
figure(fig)
hold on

plot(meshAentryWall(:,1),meshAentryWall(:,4),'blue');
plot(meshBentryWall(:,1),meshBentryWall(:,4),'red');
plot(meshCentryWall(:,1),meshCentryWall(:,4),'green');
plot(meshDentryWall(:,1),meshDentryWall(:,4),'magenta');
xlabel('x');
ylabel('Pressure, {\it p}');
legend ('Mesh A','Mesh B', 'Mesh C', 'Mesh D');
title('Spatial Distribution of Pressure at t = 4 on the Entry Wall');
hold off
%% 
mesh = ['Mesh A', 'Mesh B', 'Mesh C', 'Mesh D']'

p0 = zeros(4,1);
p0(1) = meshAentryWall(end,4);
p0(2) = meshBentryWall(end,4);
p0(3) = meshCentryWall(end,4);
p0(4) = meshDentryWall(end,4);

T = table({'Mesh A'; 'Mesh B'; 'Mesh C'; 'Mesh D'}, p0)
T.Properties.VariableNames = {'Mesh','p0'}

