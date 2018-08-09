% this is a code for modelling of concentric tube robot in free space based on " Design
% and Control of Concentric-Tube Robots " by Dupont

clearvars
clc

%% Initializing parameters 

param  % load tube parameters inside param.m file
    
tic
l=0.01*[45 30 20];   % length of tubes 
B=0.01*[-14 -10 -5];  % length of tubes before template
l_k=0.01*[10 10 15]; % length of curved part of tubes

%initial angles
alpha_1=0;
alpha_2=0;
alpha_3=pi/2;
alpha=[0 alpha_2-alpha_1 alpha_3-alpha_1];

% segmenting tubes  
% check all inputs must have n elements, n is number of tubes
[L,d_tip,EE,UUx,UUy,II,GG,JJ] = segmenting(E,Ux,Uy,I,G,J,l,B,l_k);

k=length(l); 
% figure(1)
% xlabel('S [mm]')
% hold on
% for i=1:k
%     plot(linspace(B(i),d_tip(i),10),i*ones(1,10),'r' ,'LineWidth',i*1.5)
% end
SS=L;
for i=1:length(L)
    SS(i)=sum(L(1:i));
%     plot((B(1)+SS(i))*ones(1,10),1:10,'b' ,'LineWidth',2)
end
%hold off

% S is segmented abssica of tube after template
 S=SS(SS+min(B)>0)+min(B);
 E=zeros(n,length(S)); I=E; G=E; J=E; Ux=E; Uy=E;
 for i=1:n
    E(i,:)=EE(i,SS+min(B)>0); I(i,:)=II(i,SS+min(B)>0); G(i,:)=GG(i,SS+min(B)>0);
    J(i,:)=JJ(i,SS+min(B)>0); Ux(i,:)=UUx(i,SS+min(B)>0); Uy(i,:)=UUy(i,SS+min(B)>0);
 end
 % each (i,j) element of above matrices correspond to the jth segment of
 % ith tube, 1st tube is the most inner
  
%% Solving ode for segments

span=[0 S];       % vector of tube abssica starting at zero
Length=[]; r=[];   % solved length, curvatures, and twist angles
%U1_after=[0;0;0];             % 1st tube initial curvature at segment beginning
r0=[ 0 0 0]'; R0=[cos(alpha_1) sin(alpha_1) 0; -sin(alpha_1) cos(alpha_1) 0; 0 0 1]; % why should it be R' not R?
R0=reshape(R0,[9,1]);
for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.0001];
y0_1=[r0 ; R0];

y0_2=zeros(2*n,1);
y0_2(n+1:2*n)=alpha;

y_0=[y0_2; y0_1];

[s,y] = ode45(@(s,y) ode4(s,y,Ux(:,seg),Uy(:,seg),E(:,seg),I(:,seg),G(:,seg),J(:,seg),n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i
shape=[y(:,2*n+1),y(:,2*n+2),y(:,2*n+3)];
Length=[Length; s];
r=[r; shape];

r0=shape(end,:)';
R0=y(end,2*n+4:2*n+12)';
end

[~, tube2_end] = min(abs(Length-d_tip(2)));
r2=[r(1:tube2_end,1),r(1:tube2_end,2),r(1:tube2_end,3)];
[~, tube3_end] = min(abs(Length-d_tip(3)));
r3=[r(1:tube3_end,1),r(1:tube3_end,2),r(1:tube3_end,3)];

figure(1);
plot3(r(:,1),r(:,2),r(:,3),'b','LineWidth',2)
hold on
plot3(r2(:,1),r2(:,2),r2(:,3),'r','LineWidth',4)
plot3(r3(:,1),r3(:,2),r3(:,3),'g','LineWidth',6)
xlabel('X [mm]'); ylabel('Y [m]'); zlabel('Z [m]')
grid on
axis equal

toc

