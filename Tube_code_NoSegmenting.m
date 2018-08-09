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
alpha_3=0;
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
    %plot((B(1)+SS(i))*ones(1,10),1:10,'b' ,'LineWidth',2)
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
%% Fiting smooth functions to params
span=[0 S]; 
N=50; %number of segments for linspace command
tube_length=[]; EI=zeros(n,length(S)*N); GJ=EI;
Uxx=EI; Uyy=EI;

for i=1:length(S)
    tube_length=[tube_length linspace(span(i), span(i+1)-0.000001,N)]; 
    for k=1:n
    EI(k,1+(i-1)*N:i*N)=E(k,i)*I(k,i)*ones(1,N);
    GJ(k,1+(i-1)*N:i*N)=G(k,i)*J(k,i)*ones(1,N);
    Uxx(k,1+(i-1)*N:i*N)=Ux(k,i)*ones(1,N);
    Uyy(k,1+(i-1)*N:i*N)=Uy(k,i)*ones(1,N);
    end
end
  
 %% Solving ode for segments

span=[0 S];       % vector of tube abssica starting at zero
U_x=[]; U_y=[]; 
s_span = [0 S(end)-0.0001];
y_0=zeros(2*n,1);
y_0(n+1:2*n)=alpha;

[s,y] = ode45(@(s,y) ode3(s,y,Uxx,Uyy,EI,GJ,tube_length,n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i

% calculating curvatures along x and y
ei=zeros(n,1); gj=ei; ux=ei; uy=ei;
for m=1:length(s)

% approximating physical parameters asfunction of tube length

for i=1:n
ei(i) = interp1(tube_length,EI(i,:),s(m)); % Interpolate the data set (ft,f) at time t
gj(i) = interp1(tube_length,GJ(i,:),s(m)); % Interpolate the data set (ft,f) at time t
ux(i) = interp1(tube_length,Uxx(i,:),s(m)); % Interpolate the data set (ft,f) at time t
uy(i) = interp1(tube_length,Uyy(i,:),s(m)); % Interpolate the data set (ft,f) at time t
end
Ux=ux; Uy=uy;


% calculating summation of matrices
K=zeros(3,3);SUM=zeros(3,1);
for i=1:n
    k=diag([ei(i) ei(i) gj(i)] );
    sum=[cos(y(n+i)) -sin(y(n+i)) 0; sin(y(n+i)) cos(y(n+i)) 0; 0 0 1]*k*[Ux(i); Uy(i); 0];
    K=K+k;
    SUM=SUM+sum;
end


% calculating 1st tube's curvatures in x and y direction
ux=zeros(n,1);uy=zeros(n,1);
u1= K\ SUM ;
ux(1)=u1(1); uy(1)=u1(2);

% calculating tube's curvatures in x and y direction
for i=2:n    
u= [cos(y(n+i)) sin(y(n+i)) 0; -sin(y(n+i)) cos(y(n+i)) 0; 0 0 1] * u1; 
ux(i)=u(1); uy(i)=u(2);    
end


U_x=[U_x; ux'];
U_y=[U_y; uy'];
end

Length=s;
U_z=y(:,1:n );
%Alpha=[Alpha; y(:,n+1:2*n)];




%% Calculating Shape


R0=[cos(alpha_1) -sin(alpha_1) 0; sin(alpha_1) cos(alpha_1) 0; 0 0 1]; % why should it be R' not R?

y_0=[0; 0 ;0 ;reshape(R0,[9,1])];
[s,y] = ode45(@(s,y) ode2(s,y,U_x,U_y,U_z,Length), [0 Length(end)], y_0);
%y(1) to y(3) are x , y, and z position of point materials
toc
r1=[y(:,1),y(:,2),y(:,3)];
[~, tube2_end] = min(abs(s-d_tip(2)));
r2=[y(1:tube2_end,1),y(1:tube2_end,2),y(1:tube2_end,3)];
[~, tube3_end] = min(abs(s-d_tip(3)));
r3=[y(1:tube3_end,1),y(1:tube3_end,2),y(1:tube3_end,3)];

figure(1)
xlabel('X [mm]'); ylabel('Y [m]'); zlabel('Z [m]')
plot3(r1(:,1),r1(:,2),r1(:,3),'LineWidth',2)
hold on
plot3(r2(:,1),r2(:,2),r2(:,3),'LineWidth',4)
plot3(r3(:,1),r3(:,2),r3(:,3),'LineWidth',6)
grid on
axis equal