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

[s,y] = ode23(@(s,y) ode4(s,y,Ux(:,seg),Uy(:,seg),E(:,seg),I(:,seg),G(:,seg),J(:,seg),n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i
shape=[y(:,2*n+1),y(:,2*n+2),y(:,2*n+3)];
Length=[Length; s];
r=[r; shape];

r0=shape(end,:)';
R0=y(end,2*n+4:2*n+12)';
end

r1=r;
[~, tube2_end] = min(abs(Length-d_tip(2)));
r2=[r(1:tube2_end,1),r(1:tube2_end,2),r(1:tube2_end,3)];
[~, tube3_end] = min(abs(Length-d_tip(3)));
r3=[r(1:tube3_end,1),r(1:tube3_end,2),r(1:tube3_end,3)];
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


%% System of ODEs for n tube

function dydt = ode4(s,y,Ux,Uy,E,I,G,J,n)

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
u1= K\ SUM ;
ux(1)=u1(1); uy(1)=u1(2);

% calculating tube's curvatures in x and y direction
for i=2:n    
u= [cos(y(n+i)) sin(y(n+i)) 0; -sin(y(n+i)) cos(y(n+i)) 0; 0 0 1] * u1; 
ux(i)=u(1); uy(i)=u(2);    
end

% odes for twist
for i=1:n      
    if G(i)==0
        G(i)=1; J(i)=1;  % to avoid singularity when tube doesn't exist
    end
    dydt(i)=  (  (E(i)*I(i))/(G(i)*J(i))  ) * ( ux(i)* Uy(i) -  uy(i)* Ux(i) );  % ui_z
    dydt(n+i)=  y(i)-y(1);   %alpha_i
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
