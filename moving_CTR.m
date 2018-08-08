% this is a code for modelling of concentric tube robot in free space based on " Design
% and Control of Concentric-Tube Robots " by Dupont

%% Initializing parameters 
function [r1,r2,r3] = moving_CTR(q)

param  % load tube parameters inside param.m file

% q1 o q3 are robot base movments, q3 to q6 are rbot base rotation angle.

l=0.01*[45 30 20];   % length of tubes 
l_k=0.01*[10 10 15]; % length of curved part of tubes


B=q(1:3);  % length of tubes before template
%initial angles
alpha_1=q(4);
alpha=[0 q(5)-q(4)  q(6)-q(4)];

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
% hold off

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
Length=[]; U_x=[]; U_y=[]; U_z=[]; Alpha=[];   % solved length, curvatures, and twist angles
%U1_after=[0;0;0];             % 1st tube initial curvature at segment beginning

for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.000001];
y_0=zeros(2*n,1);
y_0(n+1:2*n)=alpha;

[s,y] = ode45(@(s,y) ode(s,y,Ux(:,seg),Uy(:,seg),E(:,seg),I(:,seg),G(:,seg),J(:,seg),n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i

% calculating curvatures along x and y
for m=1:length(s)
 
K=zeros(3,3);SUM=zeros(3,1); 
for i=1:n
    k=diag([E(i,seg)*I(i,seg) E(i,seg)*I(i,seg) G(i,seg)*J(i,seg)] );
    num=[cos(y(m,n+i)) -sin(y(m,n+i)) 0; sin(y(m,n+i)) cos(y(m,n+i)) 0; 0 0 1]*k*[Ux(i,seg); Uy(i,seg); 0];
    K=K+k;
    SUM=SUM+num;
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

end


%% Calculating Shape

R0=[cos(alpha_1) sin(alpha_1) 0; -sin(alpha_1) cos(alpha_1) 0; 0 0 1]; % why should it be R' not R?

y_0=[0; 0 ;0 ;reshape(R0,[9,1])];
[s,y] = ode45(@(s,y) ode2(s,y,U_x,U_y,U_z,Length), [0 Length(end)], y_0);
%y(1) to y(3) are x , y, and z position of point materials

r1=[y(:,1),y(:,2),y(:,3)];
[~, tube2_end] = min(abs(s-d_tip(2)));
r2=[y(1:tube2_end,1),y(1:tube2_end,2),y(1:tube2_end,3)];
[~, tube3_end] = min(abs(s-d_tip(3)));
r3=[y(1:tube3_end,1),y(1:tube3_end,2),y(1:tube3_end,3)];

end
