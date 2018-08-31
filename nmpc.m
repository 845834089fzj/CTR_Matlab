
% this is a code for modle predictive control of CTR


clearvars
clc


DT=0.01;   %sampling time

r0=[0.0452679818584254,0.215646145981389,0.156739201416997]'; % initial point of robot tip--update if you change beta or alphas
r0=[[0.0131340997991548,0.109763500927478,0.152181673489298]]';

t=0:DT:10; % time interval
T_horizon=0.1;  % horizon
N=0.1/DT;  % number of samples in horizon

Mu_0=5;  % initial value of perturbation parameter

Wr=2.*eye(3,3);  % weghting function for position
Wu=0*eye(6,6);  % weghting function for input
Wuz=2*eye(3,3);  % weghting function for end curvatures

d_alpha=0.0001;   % step size for angles
d_beta=0.00025;  % step size for tube movement
d_uz=0.0001; % step size for initial curvature
%ds=0.1; % step size for slack variable

% joint limit offset
offset=0.005;
%max dist from template
max_dist=0.35;

% initial control inputs
% length of tubes before template
l=0.01*[55 30 20];   % length of tubes 
B=0.01*[-35 -15 -10];  % length of tubes before template
l_k=0.01*[10 10 15]; % length of curved part of tubes

%initial angles
alpha_1=3*pi/2;
alpha_2=pi/2;
alpha_3=pi;
% initial curvature along z
uz_init=[   1.331946135078349 1.732479299966701 0.828496580881414]';

% slack variables
s0=0.015.*[1 1 1 1]';


% param for posetive defenitness
gamma1=1; gamma2=1;
tau=0.95;
sigma=0.75;
%% estimating initial value of v and w  *note: H 4*1 dH 4*9 w 4*1 dG 1*9  dF 9*1

s=s0;
Mu=Mu_0;
rd=r0;   % desired trajectory
% calculating gradient of F and G, F=(r-r_d)^2+q^2,  G=uzi(l)^2
% no gradient
uz_0=uz_init;  % optim variable-curvature at beginning to have zero curve at ends
z=[B(1) B(2) B(3) alpha_1 alpha_2 alpha_3]; %initial value of z - control variables

[r,~,~,Uz] = moving_CTR_nmpc(z,uz_0,l,l_k); r=r(end,:)';
F0=(rd-r)'*Wr*(rd-r)+z*Wu*z';
G0=Uz'*Wuz*Uz; 
% gradient with respect to z and uz_0


delta_F= zeros(9,1); delta_G=zeros(1,9);
for p=1:9
    z_pert=z;
    uz_0_pert=uz_0;
    if p<=3
        d_z=d_beta;
        z_pert(p)=z_pert(p)+d_z;
        [r,~,~,~] = moving_CTR_nmpc(z_pert,uz_0,l,l_k); r=r(end,:)'; dF=((rd-r)'*Wr*(rd-r)+z_pert*Wu*z_pert'-F0)/d_z;delta_F(p)= dF;
        %dG=(Uz'*Wuz*Uz-G0)/d_z; delta_G(p)= dG;  
    elseif p<=6
        d_z=d_alpha;
        z_pert(p)=z_pert(p)+d_z;
        [r,~,~,~] = moving_CTR_nmpc(z_pert,uz_0,l,l_k); r=r(end,:)'; dF=((rd-r)'*Wr*(rd-r)+z_pert*Wu*z_pert'-F0)/d_z; delta_F(p)= dF;
        %dG=(Uz'*Wuz*Uz-G0)/d_z; delta_G(p)= dG;  
    else
        d_z=d_uz;
        uz_0_pert(p-6)=uz_0_pert(p-6)+d_z;
        [r,~,~,Uz] = moving_CTR_nmpc(z,uz_0_pert,l,l_k); r=r(end,:)'; %dF=((rd-r)'*Wr*(rd-r)+z*Wu*z'-F0)/d_z; delta_F(p)= dF;
        dG=(Uz'*Wuz*Uz-G0)/d_z; delta_G(p)= dG;  
    end
   
    
end


% calculating gradient of H, ..
%H= [-b3-offset>0  -b2+b3-offset>0  -b1+b2-offset>0  b1+max_dist>0];
H0=[-z(3)-offset-s(1); -z(2)+z(3)-offset-s(2); -z(1)+z(2)-offset-s(3); z(1)+max_dist-s(4)];
delta_H=[ 0 0 -1 1; 0 -1 1 0;-1 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0; 0 0 0 0]';


% computing v0 and w0
w= diag(s)\ (Mu*[1 1 1 1]');
v= (-delta_G') \ (delta_H'*w-delta_F);

%% Solving nonlinear problem


for i=1
% calculating gradient of F and G, F=(r-r_d)^2+q^2,  G=uzi(l)^2

K_max=200;
Z_optim=zeros(K_max,9);
E_optim=zeros(K_max,1);
G_optim=zeros(K_max,1);
Z_optim(1,:)=[z uz_0'];


rd=r0+[0 0 0]';   % desired trajectory
E=Mu; E2=Mu;
k=0;

while k < K_max && E2 >= Mu && sum(s>0)==4

k=k+1;
gamma=gamma1;
% d_alpha=(1/k)*0.001;   % step size for angles
% d_beta=(1/k)*0.00025;  % step size for tube movement
% d_uz=(1/k)*0.001; % step size for initial curvature


[r,~,~,Uz] = moving_CTR_nmpc(z,uz_0,l,l_k); r=r(end,:)';
F0=sqrt((rd-r)'*Wr*(rd-r))+sqrt(z*Wu*z');
G0=Uz'*Wuz*Uz; 
% gradient with respect to z and uz_0



delta_F= zeros(9,1); delta_G=zeros(1,9);
for p=1:9
    z_pert=z;
    uz_0_pert=uz_0;
    if p<=3
        d_z=d_beta;
        z_pert(p)=z_pert(p)+d_z;
        [r,~,~,Uz] = moving_CTR_nmpc(z_pert,uz_0,l,l_k); r=r(end,:)'; dF=( sqrt((rd-r)'*Wr*(rd-r))+sqrt(z_pert*Wu*z_pert')-F0)/d_z;
        delta_F(p)= dF;
        dG=(sqrt(Uz'*Wuz*Uz)-G0)/d_z; delta_G(p)= dG;  
    elseif p<=6
        d_z=d_alpha;
        z_pert(p)=z_pert(p)+d_z;
        [r,~,~,Uz] = moving_CTR_nmpc(z_pert,uz_0,l,l_k); r=r(end,:)'; dF=(sqrt((rd-r)'*Wr*(rd-r))+sqrt(z_pert*Wu*z_pert')-F0)/d_z; 
        delta_F(p)= dF;
        dG=(sqrt(Uz'*Wuz*Uz)-G0)/d_z; delta_G(p)= dG;  
    else
        d_z=d_uz;
        uz_0_pert(p-6)=uz_0_pert(p-6)+d_z;
        [r,~,~,Uz] = moving_CTR_nmpc(z,uz_0_pert,l,l_k); r=r(end,:)'; dF=(sqrt((rd-r)'*Wr*(rd-r))+sqrt(z_pert*Wu*z_pert')-F0)/d_z; 
        delta_F(p)= dF;
        dG=(sqrt(Uz'*Wuz*Uz)-G0)/d_z; delta_G(p)= dG;  
    end
   
    
end


% calculating gradient of H, ..
%H= [-b3-offset>0  -b2+b3-offset>0  -b1+b2-offset>0  b1+max_dist>0];
H0=[-z(3)-offset-s(1); -z(2)+z(3)-offset-s(2); -z(1)+z(2)-offset-s(3); z(1)+max_dist-s(4)];
H=[-z(3)-offset; -z(2)+z(3)-offset; -z(1)+z(2)-offset; z(1)+max_dist];
delta_H=[ 0 0 -1 1; 0 -1 1 0;-1 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0; 0 0 0 0]';

 delta_wH1= (w'*[-z(3)-offset-s(1); -z(2)+z(3)-offset-s(2); -z(1)-d_beta+z(2)-offset-s(3); z(1)+d_beta+max_dist-s(4)] - w'*H0 )/d_beta;
 delta_wH2= (w'*[-z(3)-offset-s(1); -z(2)-d_beta+z(3)-offset-s(2); -z(1)+z(2)+d_beta-offset-s(3); z(1)+max_dist-s(4)] - w'*H0 )/d_beta;
 delta_wH3= (w'*[-z(3)-d_beta-offset-s(1); -z(2)+z(3)+d_beta-offset-s(2); -z(1)+z(2)-offset-s(3); z(1)+max_dist-s(4)] - w'*H0 )/d_beta;
 delta_wH=[delta_wH1 delta_wH2 delta_wH3 0 0 0 0 0 0]';
%delta_wH=delta_H'*w;

% calculating hessian of L, L=F-v.G-w'*(H-s)
L0=F0-v*G0-w'*(H0-s);ddL=zeros(9,9);
for m=1:9
    for j=1:9
        ddL(m,j)= 0.5* delta_F(m) + 0.5* delta_F(j) ...
            - 0.5* v*delta_G(m) - 0.5* v* delta_G(j) - 0.5* delta_wH(m) - 0.5* delta_wH(j);
    end
end


%solveing linear system of equations for dz ds dv dw
%Sigma=Mu.*diag(s)^-2;
Sigma=(diag(s)^-1)*diag(w);

while rcond(ddL+gamma*eye(9,9))<0.01
    gamma=gamma+0.1;
end
cond=rcond(ddL+gamma*eye(9,9));

A=[ ddL+gamma*eye(9,9) zeros(9,4) -delta_G' -delta_H'; ...
    zeros(4,9) Sigma zeros(4,1) eye(4,4) ;...
   -delta_G   zeros(1,4) gamma2  zeros(1,4); ...
   -delta_H eye(4,4) zeros(4,1) zeros(4,4)];
b=-[ delta_F-delta_G'*v-delta_H'*w  ; w- Mu*(diag(s)^-1)*[1 1 1 1]'; -G0; -H+s];

d=A\b;

dz=d(1:9); ds=d(10:13); dv=d(14); dw=d(15:18);

% stopping criteria
E2=max(b.^2)
E=F0+G0
% find new iteration

% finding new step sizes
cost=-1;a_w=1;
while(sum(cost<0)>0)
cost= tau*w+a_w*dw ;
a_w=a_w*0.95;
end
cost=-1;a_s=1;
while(sum(cost<0)>0)
cost= tau*s+a_s*ds ;
a_s=a_s*0.95;
end


z=z+a_s.*dz(1:6)';
uz_0=uz_0+a_s.*dz(7:9);
s=s+a_s.*ds; %s(s<0)=0;
v=v+a_w.*dv;
w=w+a_w.*dw;

Mu= sigma*Mu;

E_optim(k)=E;
G_optim(k)=G0;     
Z_optim(k+1,:)=[z uz_0'];


end

[~,index]=min(E_optim(E_optim>0));
E_optim(index)-G_optim(index)
G_optim(index)
z=Z_optim(index,1:6);
uz_0=Z_optim(index,7:9)';
k


[r,r2,r3,Uz] = moving_CTR_nmpc(z,uz_0);

figure(1);
plot3(r(:,1),r(:,2),r(:,3),'k','LineWidth',2)
hold on
plot3(r2(:,1),r2(:,2),r2(:,3),'k','LineWidth',4)
plot3(r3(:,1),r3(:,2),r3(:,3),'k','LineWidth',6)
xlabel('X [mm]'); ylabel('Y [m]'); zlabel('Z [m]')
grid on
axis equal

sqrt((rd-r(end,:)')'*Wr*(rd-r(end,:)') )
sqrt(Uz'*Wuz*Uz)

end


%%

function [r1,r2,r3,Uz] = moving_CTR_nmpc(q,uz_0,l,l_k)

param  % load tube parameters inside param.m file

% q1 o q3 are robot base movments, q3 to q6 are rbot base rotation angle.


B=q(1:3);  % length of tubes before template
%initial angles
alpha_1=q(4);
alpha=[q(4) q(5) q(6)];


% segmenting tubes  
% check all inputs must have n elements, n is number of tubes
[L,d_tip,EE,UUx,UUy,II,GG,JJ] = segmenting(E,Ux,Uy,I,G,J,l,B,l_k);

SS=L;
for i=1:length(L)
    SS(i)=sum(L(1:i));
%     plot((B(1)+SS(i))*ones(1,10),1:10,'b' ,'LineWidth',2)
end

% S is segmented abssica of tube after template
 S=SS(SS+min(B)>0)+min(B);
 E=zeros(n,length(S)); I=E; G=E; J=E; Ux=E; Uy=E;
 for i=1:n
    E(i,:)=EE(i,SS+min(B)>0); I(i,:)=II(i,SS+min(B)>0); G(i,:)=GG(i,SS+min(B)>0);
    J(i,:)=JJ(i,SS+min(B)>0); Ux(i,:)=UUx(i,SS+min(B)>0); Uy(i,:)=UUy(i,SS+min(B)>0);
 end
 % each (i,j) element of above matrices correspond to the jth segment of
 % ith tube, 1st tube is the most inner

 %% Solving ode for shape

span=[0 S];       % vector of tube abssica starting at zero
Length=[]; r=[]; U_z=[]; angle=[]; % solved length, curvatures, and twist angles
%U1_after=[0;0;0];             % 1st tube initial curvature at segment beginning
r0=[ 0 0 0]'; R0=[cos(alpha_1) sin(alpha_1) 0; -sin(alpha_1) cos(alpha_1) 0; 0 0 1];
R0=reshape(R0,[9,1]);
alpha=alpha-B.*uz_0'; 

for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.0000001];
y0_1=[r0 ; R0];

y0_2=zeros(2*n,1);
y0_2(n+1:2*n)=alpha;
y0_2(1:n)=uz_0;

y_0=[y0_2; y0_1];

[s,y] = ode45(@(s,y) ode5(s,y,Ux(:,seg),Uy(:,seg),E(:,seg),I(:,seg),G(:,seg),J(:,seg),n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i
shape=[y(:,2*n+1),y(:,2*n+2),y(:,2*n+3)];
Length=[Length; s];
r=[r; shape];
U_z=[U_z; y(:,1:n )];
r0=shape(end,:)';
R0=y(end,2*n+4:2*n+12)';
angle=[angle; y(:,1+n:2*n )];
% 
 uz_0=U_z(end,:)';
 alpha=[y(end,n+1),y(end,n+2),y(end,n+3)]';

end

Uz=zeros(n,1);
for i=1:n
[~,index] =  min( abs(Length-d_tip(i)+0.0001) );
Uz(i)= U_z(index,i);
end

r1=r;
[~, tube2_end] = min(abs(Length-d_tip(2)));
r2=[r(1:tube2_end,1),r(1:tube2_end,2),r(1:tube2_end,3)];
[~, tube3_end] = min(abs(Length-d_tip(3)));
r3=[r(1:tube3_end,1),r(1:tube3_end,2),r(1:tube3_end,3)];

end



%% ODE
function dydt = ode5(~,y,Ux,Uy,E,I,G,J,n)

dydt=zeros(2*n+12,1);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% second n elements of y are twist angles, alpha_i
% last 12 elements are r (position) and R (orientations), respectively

% calculating summation of matrices
K=zeros(3,3);SUM=zeros(3,1);
for i=1:n
    k=diag([E(i)*I(i) E(i)*I(i) G(i)*J(i)] );
    sum=[cos(y(n+i)) -sin(y(n+i)) 0; sin(y(n+i)) cos(y(n+i)) 0; 0 0 1]*k*[Ux(i); Uy(i); 0];
    K=K+k;
    SUM=SUM+sum;
end


% calculating 1st tube's curvatures in x and y direction
ux=zeros(n,1);uy=zeros(n,1);

% calculating tube's curvatures in x and y direction
for i=1:n    
u= K\([cos(y(n+i)) sin(y(n+i)) 0; -sin(y(n+i)) cos(y(n+i)) 0; 0 0 1] * SUM ); 
ux(i)=u(1); uy(i)=u(2);    
end

% odes for twist
for i=1:n      
    if G(i)==0
        G(i)=1; J(i)=1;  % to avoid singularity when tube doesn't exist
    end
    dydt(i)=  (  (E(i)*I(i))/(G(i)*J(i))  ) * ( ux(i)* Uy(i) -  uy(i)* Ux(i) );  % ui_z
    dydt(n+i)=  y(i);   %alpha_i
end


e3=[0 0 1]';              
uz = y(1:n); 



% y(1) to y(3) are position of point materials
%r1=[y(1); y(2); y(3)];
% y(4) to y(12) are rotation matrix elements
R1=[y(2*n+4) y(2*n+5) y(2*n+6);y(2*n+7) y(2*n+8) y(2*n+9);y(2*n+10) y(2*n+11) y(2*n+12)];


u_hat=[0 -uz(1) uy(1) ; uz(1) 0 -ux(1) ; -uy(1) ux(1) 0 ];


% odes
dr1 = R1*e3;
dR1=R1*u_hat;


dydt(2*n+1)=dr1(1);dydt(2*n+2)=dr1(2);dydt(2*n+3)=dr1(3);
dR=dR1'; 
dR=dR(:);
for i=4:12
   dydt(2*n+i)=dR(i-3);
end

end



%% code for segmenting tubes

function [L,d1,E,Ux,Uy,I,G,J] = segmenting(E,Ux,Uy,I,G,J,l,B,l_k)

% all vectors must be sorted, starting element belongs to the most inner tube
%E, U, I, G, J   stifness, curvature, inertia, torsion constant, and second moment of inertia vectors for each tube
%l vector of tube length
%B  vector of tube movments with respect to template position, i.e., s=0 (always negative)
%l_k vecot oftube's curved part length

k=length(l); 

d1= l+B; % position of tip of the tubes
d2=d1-l_k; % position of the point where tube bending starts
points=[0 B d2 d1];
[L, index]=sort(points);
L = 1e-5*floor(1e5*diff(L));  % length of each segment 
%(used floor because diff command doesn't give absolute zero sometimes)

for i=1:k-1
if B(i)>B(i+1)
    sprintf('inner tube is clashing into outer tubes')
    E=zeros(k,length(L));
    I=E; G=E; J=E; Ux=E; Uy=E;
    return
end
end

EE=zeros(k,length(L));
II=EE; GG=EE; JJ=EE; UUx=EE; UUy=EE;

for i=1:k
    
a=find(index==i+1);   % find where tube begins
b=find(index==1*k+i+1); % find where tube curve starts
c=find(index==2*k+i+1); % find where tube ends

if L(a)==0; a=a+1;  end
if L(b)==0; b=b+1;  end
if c<=length(L)
    if L(c)==0; c=c+1; end
end
    
EE(i,a:c-1)=E(i);
II(i,a:c-1)=I(i);
GG(i,a:c-1)=G(i);
JJ(i,a:c-1)=J(i);
UUx(i,b:c-1)=Ux(i);
UUy(i,b:c-1)=Uy(i);
end

l=L(~(L==0));  % get rid of zero lengthes
E=zeros(k,length(l)); I=E; G=E; J=E; Ux=E; Uy=E;
 for i=1:k
    E(i,:)=EE(i,~(L==0)); I(i,:)=II(i,~(L==0)); G(i,:)=GG(i,~(L==0));
    J(i,:)=JJ(i,~(L==0)); Ux(i,:)=UUx(i,~(L==0)); Uy(i,:)=UUy(i,~(L==0));
 end
L=L(~(L==0));

end


