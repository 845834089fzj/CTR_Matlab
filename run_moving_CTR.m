
clearvars
clc

t=0:0.01:10;
v1=0.5; v2=0.3; v3=0.2;
w1=0.05; w2=0.03; w3=0.05;

for i=1:length(t)
% length of tubes before template
B=0.01*[-14+v1*t(i) -10+v2*t(i) -5+v3*t(i)]; 
%initial angles
alpha_1=3*pi/2+w1*t(i);
alpha_2=pi+w1*t(i);
alpha_3=pi+w1*t(i);
q=[B alpha_1 alpha_2 alpha_3];
    
tic
[r1,r2,r3] = moving_CTR(q);
toc
% figure(2)
% clf(figure(2))
% hold on
% plot3(r1(:,1),r1(:,2),r1(:,3),'b','LineWidth',2)
% plot3(r2(:,1),r2(:,2),r2(:,3),'r','LineWidth',4)
% plot3(r3(:,1),r3(:,2),r3(:,3),'g','LineWidth',6)
% xlabel('X [mm]'); ylabel('Y [m]'); zlabel('Z [m]')
% grid on
% axis equal

end