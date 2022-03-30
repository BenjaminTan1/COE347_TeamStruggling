%Vertex creator
r=.5;%radius of cylinder
r2=1;%radius of boundary layer
upper_height=8*r;%upper height=lower height
forward_length=-4;
backward_length=12;

%establish vertice locations
vertices=zeros(64,3);
%inner circle
for z=1:2
    for i=1:8
        theta=(i-1)*pi/4;
        vertices(i+32*(z-1),:)=[r*cos(theta),r*sin(theta),(z-1)*.1];
        vertices(i+8+32*(z-1),:)=[r2*cos(theta),r2*sin(theta),(z-1)*.1];
        
    end
    up1=r2*sin(pi/4);
    down1=-up1;
    right1=r2*cos(pi/4);
    left1=-right1;
    vertices(17+32*(z-1),:)=[backward_length,0,(z-1)*.1];
    vertices(18+32*(z-1),:)=[backward_length,up1,(z-1)*.1];
    vertices(19+32*(z-1),:)=[backward_length,upper_height,(z-1)*.1];
    vertices(20+32*(z-1),:)=[right1,upper_height,(z-1)*.1];
    vertices(21+32*(z-1),:)=[0,upper_height,(z-1)*.1];
    vertices(22+32*(z-1),:)=[left1,upper_height,(z-1)*.1];
    vertices(23+32*(z-1),:)=[forward_length,upper_height,(z-1)*.1];
    vertices(24+32*(z-1),:)=[forward_length,up1,(z-1)*.1];
    vertices(25+32*(z-1),:)=[forward_length,0,(z-1)*.1];
    vertices(26+32*(z-1),:)=[forward_length,down1,(z-1)*.1];
    vertices(27+32*(z-1),:)=[forward_length,-upper_height,(z-1)*.1];
    vertices(28+32*(z-1),:)=[left1,-upper_height,(z-1)*.1];
    vertices(29+32*(z-1),:)=[0,-upper_height,(z-1)*.1];
    vertices(30+32*(z-1),:)=[right1,-upper_height,(z-1)*.1];
    vertices(31+32*(z-1),:)=[backward_length,-upper_height,(z-1)*.1];
    vertices(32+32*(z-1),:)=[backward_length,down1,(z-1)*.1];
    
end

%print out stuff and save it in txt file
fileID=fopen('verticesA.txt','w');
for a=1:64
    fprintf(fileID,'(%e %e %e) \n',vertices(a,:));
end
fclose(fileID);
fileID=fopen('edgesA.txt','w');
%generate edges
for a =1:8
    theta=pi/4*(a-1)+pi/8;
    xP1=r*cos(theta);
    yP1=r*sin(theta);
    xP2=r2*cos(theta);
    yP2=r2*sin(theta);
    z1=0;
    z2=.1;
    fprintf(fileID,'arc %d %d (%e %e %e) \n',[a-1,a,xP1, yP1,z1]);
    fprintf(fileID,'arc %d %d (%e %e %e) \n', [a+7,a+8,xP2,yP2,z1]);
    fprintf(fileID,'arc %d %d (%e %e %e) \n', [a+31,a+32,xP1,yP1,z2]);
    fprintf(fileID, 'arc %d %d (%e %e %e) \n', [a+39, a+40, xP2, yP2,z2]);
    
    


end

figure(1)

scatter(vertices(:,1),vertices(:,2), 'filled')
hold on
for i=1:32

text(vertices(i,1)+.02,vertices(i,2)+.05,num2str(i-1),'FontSize',8)
text(vertices(i+32,1)+.175,vertices(i+32,2)+.05,num2str(i+31),'FontSize',8, 'Color',[1,0,0])
end
title('Vertex Schematic B')
xlabel('x (m)')
ylabel('y (m)')