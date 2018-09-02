% this is a code for modelling of concentric tube robot in free space based on " Design
% and Control of Concentric-Tube Robots " by Dupont

clearvars
clc

%% Initializing parameters 
global Ux Uy E I G J n  S  span a d_tip B

param  % load tube parameters inside param.m file
    

l=0.01*[55 30 20];   % length of tubes 
B=0.01*[-35 -15 -10];  % length of tubes before template
l_k=0.01*[10 10 15]; % length of curved part of tubes

%initial angles % direction of rotation for 1st tube seems to be different
%from 2nd and 3rd
alpha_1=pi/2;
alpha_2=pi;
alpha_3=pi/3;
a=[0 alpha_2-alpha_1 alpha_3-alpha_1];
    
Ux=[-5 5 5];
Uy=[10 8 2];

% segmenting tubes  
% check all inputs must have n elements, n is number of tubes
[L,d_tip,EE,UUx,UUy,II,GG,JJ] = segmenting(E,Ux,Uy,I,G,J,l,B,l_k);

k=length(l); 
figure(1)
xlabel('S [mm]')
hold on
for i=1:k
    plot(linspace(B(i),d_tip(i),10),i*ones(1,10),'r' ,'LineWidth',i*1.5)
end
SS=L;
for i=1:length(L)
    SS(i)=sum(L(1:i));
    plot((B(1)+SS(i))*ones(1,10),1:10,'b' ,'LineWidth',2)
end
hold off

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
% span=[0 S]; 
% N=50; %number of segments for linspace command
% tube_length=[]; EI=zeros(n,length(S)*N); GJ=EI;
% Uxx=EI; Uyy=EI;
% 
% for i=1:length(S)
%     tube_length=[tube_length linspace(span(i), span(i+1)+i*0.001,N)]; 
%     for k=1:n
%     EI(k,1+(i-1)*N:i*N)=E(k,i)*I(k,i)*ones(1,N);
%     GJ(k,1+(i-1)*N:i*N)=G(k,i)*J(k,i)*ones(1,N);
%     Uxx(k,1+(i-1)*N:i*N)=Ux(k,i)*ones(1,N);
%     Uyy(k,1+(i-1)*N:i*N)=Uy(k,i)*ones(1,N);
%     end
% end
%   
 %% Solving ode for segments

span=[0 S];       % vector of tube abssica starting at zero
%U1_after=[0;0;0];             % 1st tube initial curvature at segment beginning

x0=[0 0 0]'; A=[]; b=[];
x = fmincon(@Cost,x0,A,b)

y_0=zeros(2*n,1);
y_0(n+1:2*n)=a;
y_0(1:n)=x;
[U_x, U_y, U_z, Alpha, Length, C] = curvature(Ux,Uy,E,I,G,J,n, S, span,y_0);
%% Calculating Shape

%R0=eye(3,3);
R0=[cos(alpha_1) sin(alpha_1) 0; -sin(alpha_1) cos(alpha_1) 0; 0 0 1];

y_0=[0; 0 ;0 ;reshape(R0,[9,1])];
[s,y] = ode45(@(s,y) ode2(s,y,U_x,U_y,U_z,Length), [0 Length(end)], y_0);
%y(1) to y(3) are x , y, and z position of point materials

figure(2)
plot3(y(:,1),y(:,2),y(:,3), 'LineWidth',2)
xlabel('X [mm]'); ylabel('Y [m]'); zlabel('Z [m]')
grid on
axis equal
hold on
%% Solver for finding cost function (value of Uz at tip)


function C= Cost(uz_0)
global Ux Uy E I G J n  S span a d_tip B

y_0=zeros(2*n,1);
y_0(n+1:2*n)=a'-B'.*uz_0;
y_0(1:n)=uz_0;

C=0;
U_z=[]; 
for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.00001];


[s,y] = ode45(@(s,y) ode(s,y,Ux(:,seg),Uy(:,seg),E(:,seg),I(:,seg),G(:,seg),J(:,seg),n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i

U_z=[U_z; y(:,1:n )];

% check for BC at end of tubes, U_z should be zero
for p=1:n
if (abs(s(end)-d_tip(p))<0.001)
C=1e8.*U_z(end,:)*U_z(end,:)'+C;
end
end

end
C
end


%% Solver for tube curvature

function [U_x, U_y, U_z, Alpha, Length, C] = curvature(Ux,Uy,E,I,G,J,n, S, span,y_0)

Length=[]; U_x=[]; U_y=[]; U_z=[]; Alpha=[];   % solved length, curvatures, and twist angles
for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.0001];


[s,y] = ode45(@(s,y) ode(s,y,Ux(:,seg),Uy(:,seg),E(:,seg),I(:,seg),G(:,seg),J(:,seg),n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i

% calculating curvatures along x and y
for m=1:length(s)
 
K=zeros(3,3);SUM=zeros(3,1);
for i=1:n
    k=diag([E(i,seg)*I(i,seg) E(i,seg)*I(i,seg) G(i,seg)*J(i,seg)] );
    sum=[cos(y(m,n+i)) -sin(y(m,n+i)) 0; sin(y(m,n+i)) cos(y(m,n+i)) 0; 0 0 1]*k*[Ux(i,seg); Uy(i,seg); 0];
    K=K+k;
    SUM=SUM+sum;
end

ux=zeros(1,n);uy=zeros(1,n);
u1= K\ SUM;
ux(1)=u1(1); uy(1)=u1(2);

for i=2:n    
u= [cos(y(m,n+i)) sin(y(m,n+i)) 0; -sin(y(m,n+i)) cos(y(m,n+i)) 0; 0 0 1] * u1;
ux(i)=u(1); uy(i)=u(2);    
end

U_x=[U_x; ux];
U_y=[U_y; uy];
end

Length=[Length; s];
U_z=[U_z; y(:,1:n )];
Alpha=[Alpha; y(:,n+1:2*n)];

% calculating boundary conditions
% m_before=0; m_after=0;
% for i=1:n
%     k_before=diag([E(i,seg)*I(i,seg) E(i,seg)*I(i,seg) G(i,seg)*J(i,seg)] ); 
%     k_after=diag([E(i,seg+1)*I(i,seg+1) E(i,seg+1)*I(i,seg+1) G(i,seg+1)*J(i,seg+1)] );
%     U_before=[U_x(end,i); U_y(end,i); U_z(end,i)];
%     U_star_before=[ Ux(i,seg); Uy(i,seg); 0 ];
%     U_star_after=[ Ux(i,seg+1); Uy(i,seg+1); 0 ];
%     R_alpha=[cos(Alpha(end,i)) -sin(Alpha(end,i)) 0; sin(Alpha(end,i)) cos(Alpha(end,i)) 0; 0 0 1];
%     m_before=m_before+R_alpha*k_before* ( U_before-U_star_before );
%     m_after=m_after-R_alpha*k_after*U_star_after;
% end
% 
% U1_after= K\(m_before-m_after);


end

C=U_z(end,:)*U_z(end,:)';
end