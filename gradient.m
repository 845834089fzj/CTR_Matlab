
clearvars
clc

l=0.01*[55 30 20];   % length of tubes 
B=0.01*[-35 -15 -10];  % length of tubes before template
l_k=0.01*[10 10 15]; % length of curved part of tubes

%initial angles
alpha_1=3*pi/2;
alpha_2=pi/2;
alpha_3=pi;

uz_init=[0 0 0];

uz_0=uz_init; 
z=[B(1) B(2) B(3) alpha_1 alpha_2 alpha_3];



[r,r2,r3,Uz,a,uuz,U_x,U_y] =  moving_CTR_nmpc(z,uz_0',l,l_k); 

m=3;
z(m)=z(m)+0.0001;
%uz_0(m-6)=uz_0(m-6)+0.0001;
[dr,~,~,dUz,dda,duuz,dU_x,dU_y] =  moving_CTR_nmpc(z,uz_0',l,l_k); 



(duuz(end,:)-uuz(end,:))./0.0001
(dda(end,:)-a(end,:))./0.0001;
(dr(1,:)-r(1,:))./0.0001;


[dr,~,~,~,Dr, Duz, da,duz] = CTR_gradient(z,uz_0,l,l_k); 

ddUx=(dU_x(end,:)-U_x(end,:))./0.0001;
ddUy=(dU_x(end,:)-U_x(end,:))./0.0001;
(duuz(end,:)-uuz(end,:))./0.0001;

da(end,m:9:end);
duz(end,m:9:end)

Dr(5,:);

figure(1);
hold on
plot3(r(:,1),r(:,2),r(:,3),'k','LineWidth',2)
plot3(dr(:,1),dr(:,2),dr(:,3),'r','LineWidth',2)
grid on
axis equal

% [~,~,~,~,Dr, Duz] = CTR_gradient(z,uz_0,l,l_k); 
% 
% grad=zeros(9,1);
% for i=1:9
%     dUz=Duz(i,:)';
%     dr=Dr(i,:)';
%     grad(i)=-2*dr'*Wr*(rd-r)-2*dUz'*Wr*Uz;
% end


%%

function [r1,r2,r3,Uz, Dr, Duz,dalpha,duz] = CTR_gradient(q,uz_0,l,l_k)

param  % load tube parameters inside param.m file

% q1 o q3 are robot base movments, q3 to q6 are rbot base rotation angle.

B=q(1:3);  % length of tubes before template
%initial angles
alpha_1=q(4);
alpha=[q(4) q(5) q(6)];

%d_tip= l+B; % position of tip of the tubes


 %% Solving ode for shape
%span=linspace(0,l(1)+B(1),100);       % vector of tube abssica starting at zero

span=[0 l(1)+B(1)];       % vector of tube abssica starting at zero
r0=[ 0 0 0]'; R0=[cos(alpha_1) -sin(alpha_1) 0; sin(alpha_1) cos(alpha_1) 0; 0 0 1];
R0=reshape(R0',[9,1]);
%alpha=alpha-B.*uz_0'; 
dR0=[-sin(alpha_1) -cos(alpha_1) 0; cos(alpha_1) -sin(alpha_1) 0; 0 0 0];
dR0=reshape(dR0',[9,1]);


y_0=zeros(2*n+12+9*n+9*n+9*12,1);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% second n elements of y are twist angles, alpha_i
% next 12 elements are r (position) and R (orientations), respectively
% next 9*n elements of y are derivatives of gradient of curvatures along z 
% with respect to inputs and uz0, e.g., y= [ d(du1_z/dB(1))/ds ... ]
% next 9*12 elements are r (position) and R (orientations) gradients, respectively

y0_1=[r0 ; R0];
y0_2=zeros(2*n,1);y0_2(n+1:2*n)=alpha; y0_2(1:n)=uz_0;
duz1_0=zeros(9,1);duz1_0(7)=1; duz2_0=zeros(9,1);duz2_0(8)=1;
duz3_0=zeros(9,1);duz3_0(9)=1;
dalpha1_0=zeros(9,1);dalpha1_0(4)=1; dalpha2_0=zeros(9,1);dalpha2_0(5)=1;
dalpha3_0=zeros(9,1);dalpha3_0(6)=1;

y_0(1:2*n)=y0_2;   %uz_0 and alpha0 
y_0(2*n+1:2*n+12)=y0_1; %r0 and R0
y_0(2*n+12+1:2*n+12+n*9)=[duz1_0; duz2_0; duz3_0];
y_0(2*n+12+n*9+1:2*n+12+n*9+n*9)=[dalpha1_0; dalpha2_0; dalpha3_0];
y_0(2*n+12+2*n*9+9*3+1:2*n+12+2*n*9+9*3+9)=dR0;


[s,y] = ode45(@(s,y) gradient_ode(s,y,Ux,Uy,E.*I,G.*J,n,l,l_k,B), span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i
r=[y(:,2*n+1),y(:,2*n+2),y(:,2*n+3)];
U_z=[y(:,1:n )];
Duz=zeros(9,3); % rows are inputs z, column are tube nums
Uz=zeros(n,1);
for i=1:n
[~,index] =  min( abs(s-l(i)-B(i)) );
Uz(i)= U_z(index,i);
Duz(:,i)=y(index,2*n+12+9*(i-1)+1:2*n+12+9*i);
end

Dr=zeros(9,3); % rows are inputs z, column are x, y, and z
Dr(:,1)=y(end,2*n+12+9*n+9*n+9*9+1:3:end);
Dr(:,2)=y(end,2*n+12+9*n+9*n+9*9+2:3:end);
Dr(:,3)=y(end,2*n+12+9*n+9*n+9*9+3:3:end);

dalpha=y(:,2*n+12+9*n+1:2*n+12+9*n+9*n);
duz=y(:,2*n+12+1:2*n+12+9*n);
r1=r;
[~, tube2_end] = min(abs(s-l(2)-B(2)));
r2=[r(1:tube2_end,1),r(1:tube2_end,2),r(1:tube2_end,3)];
[~, tube3_end] = min(abs(s-l(3)-B(3)));
r3=[r(1:tube3_end,1),r(1:tube3_end,2),r(1:tube3_end,3)];

end


%% ODE
function dydt = gradient_ode(s,y,Ux,Uy,EI,GJ,n,l,l_k,B)

dydt=zeros(2*n+12+9*n+9*n+9*12,1);

% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% second n elements of y are twist angles, alpha_i
% last 12 elements are r (position) and R (orientations), respectively
% next 9*n elements of y are derivatives of gradient of curvatures along z 
% with respect to inputs and uz0, e.g., y= [ d(du1_z/dB(1))/ds ... ]
% next 9*n elements of y are derivatives of gradient of alpha
% next 9*12 elements are r (position) and R (orientations) gradients, respectively
eps=0.0005; dUx=zeros(3,1); dUy=zeros(3,1);
for i=1:n
Ux(i)=Ux(i)*0.5*(1+tanh((s- l(i)-B(i)+l_k(i))/eps));
Uy(i)=Uy(i)*0.5*(1+tanh((s- l(i)-B(i)+l_k(i))/eps));
EI(i)=EI(i)*heaviside(-s+ l(i)+B(i));

%derivative of curvatures with repect to inputs (beta)

dUx(i)=Ux(i)*dirac(s-l(i)-B(i)+l_k(i));
dUy(i)=Uy(i)*dirac(s-l(i)-B(i)+l_k(i));
dEI(i)=EI(i)*dirac(-s+ l(i)+B(i)); %*heaviside(-s+ l(i)+B(i));
dGJ(i)=GJ(i)*dirac(-s+ l(i)+B(i)); %*heaviside(-s+ l(i)+B(i));

end
% calculating 1st tube's curvatures in x and y direction

ux=zeros(n,1);
uy=zeros(n,1);
% calculating tube's curvatures in x and y direction
% calculating tube's curvatures in x and y direction
for i=1:n  
ux(i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      EI(1)*Ux(1)*cos(y(n+i)-y(n+1))+ EI(1)*Uy(1)*sin(y(n+i)-y(n+1))  + ...
      EI(2)*Ux(2)*cos(y(n+i)-y(n+2))+ EI(2)*Uy(2)*sin(y(n+i)-y(n+2))  + ...
      EI(3)*Ux(3)*cos(y(n+i)-y(n+3))+ EI(3)*Uy(3)*sin(y(n+i)-y(n+3)) );
uy(i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      -EI(1)*Ux(1)*sin(y(n+i)-y(n+1))+ EI(1)*Uy(1)*cos(y(n+i)-y(n+1))  + ...
      -EI(2)*Ux(2)*sin(y(n+i)-y(n+2))+ EI(2)*Uy(2)*cos(y(n+i)-y(n+2))  + ...
      -EI(3)*Ux(3)*sin(y(n+i)-y(n+3))+ EI(3)*Uy(3)*cos(y(n+i)-y(n+3)) ); 
end


  
% dux(j,i) is derivative of gradient of tube i curvature with respect to input j
dux=zeros(9,3);
duy=zeros(9,3);
duz=zeros(9,3);

% derivative with respect to betas (tubes translation) j=1:3
for j=1:3  % betas
    for i=1:3 % tube number
        dux(j,i)= (1/(EI(1)+EI(2)+EI(3)))* ...
            (EI(j)*dUx(j)*cos(y(n+i)-y(n+j))+ EI(j)*dUy(j)*sin(y(n+i)-y(n+j))) +...
            (-dEI(j)/(EI(1)+EI(2)+EI(3))^2)* (...
        EI(1)*Ux(1)*cos(y(n+i)-y(n+1))+ EI(1)*Uy(1)*sin(y(n+i)-y(n+1))  + ...
        EI(2)*Ux(2)*cos(y(n+i)-y(n+2))+ EI(2)*Uy(2)*sin(y(n+i)-y(n+2))  + ...
        EI(3)*Ux(3)*cos(y(n+i)-y(n+3))+ EI(3)*Uy(3)*sin(y(n+i)-y(n+3)) ) + ...
        (1/(EI(1)+EI(2)+EI(3)))* (...
        dEI(j)*Ux(j)*cos(y(n+i)-y(n+j))+ dEI(j)*Uy(j)*sin(y(n+i)-y(n+j)) )+ ...
        (1/(EI(1)+EI(2)+EI(3)))* (...
     - EI(1)*Ux(1)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+1))+ ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+1))  + ...
       EI(1)*Ux(1)*y(2*n+12+9*n+(1-1)*9+j)*sin(y(n+i)-y(n+1)) - ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(1-1)*9+j)*cos(y(n+i)-y(n+1))  + ...
       -EI(2)*Ux(2)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+2))+ ...
       EI(2)*Uy(2)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+2))  + ...
       EI(2)*Ux(2)*y(2*n+12+9*n+(2-1)*9+j)*sin(y(n+i)-y(n+2))- ...
       EI(2)*Uy(2)*y(2*n+12+9*n+(2-1)*9+j)*cos(y(n+i)-y(n+2))  + ...
       -EI(3)*Ux(3)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+3))+ ...
       EI(3)*Uy(3)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+3)) + ...
       EI(3)*Ux(3)*y(2*n+12+9*n+(3-1)*9+j)*sin(y(n+i)-y(n+3))- ...
       EI(3)*Uy(3)*y(2*n+12+9*n+(3-1)*9+j)*cos(y(n+i)-y(n+3)) );
 
        duy(j,i)= (1/(EI(1)+EI(2)+EI(3)))* ...
       (-EI(j)*dUx(j)*sin(y(n+i)-y(n+j))+ EI(j)*dUy(j)*cos(y(n+i)-y(n+j)))+...
       (-dEI(j)/(EI(1)+EI(2)+EI(3))^2)* (...
      -EI(1)*Ux(1)*sin(y(n+i)-y(n+1))+ EI(1)*Uy(1)*cos(y(n+i)-y(n+1))  + ...
      -EI(2)*Ux(2)*sin(y(n+i)-y(n+2))+ EI(2)*Uy(2)*cos(y(n+i)-y(n+2))  + ...
      -EI(3)*Ux(3)*sin(y(n+i)-y(n+3))+ EI(3)*Uy(3)*cos(y(n+i)-y(n+3)) )+...
      (1/(EI(1)+EI(2)+EI(3)))* (...
      -dEI(j)*Ux(j)*sin(y(n+i)-y(n+j))+ dEI(j)*Uy(j)*cos(y(n+i)-y(n+j)) )+...
      (1/(EI(1)+EI(2)+EI(3)))* (...
      -EI(1)*Ux(1)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+1))+ ...
      -EI(1)*Uy(1)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+1))  + ...
       EI(1)*Ux(1)*y(2*n+12+9*n+(1-1)*9+j)*cos(y(n+i)-y(n+1))+ ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(1-1)*9+j)*sin(y(n+i)-y(n+1))  + ...
      -EI(2)*Ux(2)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+2))+ ...
      -EI(2)*Uy(2)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+2))  + ...
      EI(2)*Ux(2)*y(2*n+12+9*n+(2-1)*9+j)*cos(y(n+i)-y(n+2))+ ...
      EI(2)*Uy(2)*y(2*n+12+9*n+(2-1)*9+j)*sin(y(n+i)-y(n+2))  + ...
      -EI(3)*Ux(3)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+3))+ ...
      -EI(3)*Uy(3)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+3)) + ...
      EI(3)*Ux(3)*y(2*n+12+9*n+(3-1)*9+j)*cos(y(n+i)-y(n+3))+ ...
      EI(3)*Uy(3)*y(2*n+12+9*n+(3-1)*9+j)*sin(y(n+i)-y(n+3)) ); 
    end  
end

dux(1,1);
% derivative with respect to alphas
for j=4:6  % alphas
    for i=1:3 % tube number
dux(j,i)= (1/(EI(1)+EI(2)+EI(3)))* (...
     - EI(1)*Ux(1)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+1))+ ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+1))  + ...
       EI(1)*Ux(1)*y(2*n+12+9*n+(1-1)*9+j)*sin(y(n+i)-y(n+1))- ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(1-1)*9+j)*cos(y(n+i)-y(n+1))  + ...
       -EI(2)*Ux(2)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+2))+ ...
       EI(2)*Uy(2)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+2))  + ...
       EI(2)*Ux(2)*y(2*n+12+9*n+(2-1)*9+j)*sin(y(n+i)-y(n+2))- ...
       EI(2)*Uy(2)*y(2*n+12+9*n+(2-1)*9+j)*cos(y(n+i)-y(n+2))  + ...
       -EI(3)*Ux(3)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+3))+ ...
       EI(3)*Uy(3)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+3)) + ...
       EI(3)*Ux(3)*y(2*n+12+9*n+(3-1)*9+j)*sin(y(n+i)-y(n+3))- ...
       EI(3)*Uy(3)*y(2*n+12+9*n+(3-1)*9+j)*cos(y(n+i)-y(n+3)) );
   
   duy(j,i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      -EI(1)*Ux(1)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+1))+ ...
      -EI(1)*Uy(1)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+1))  + ...
       EI(1)*Ux(1)*y(2*n+12+9*n+(1-1)*9+j)*cos(y(n+i)-y(n+1))+ ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(1-1)*9+j)*sin(y(n+i)-y(n+1))  + ...
      -EI(2)*Ux(2)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+2))+ ...
      -EI(2)*Uy(2)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+2))  + ...
      EI(2)*Ux(2)*y(2*n+12+9*n+(2-1)*9+j)*cos(y(n+i)-y(n+2))+ ...
      EI(2)*Uy(2)*y(2*n+12+9*n+(2-1)*9+j)*sin(y(n+i)-y(n+2))  + ...
      -EI(3)*Ux(3)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+3))+ ...
      -EI(3)*Uy(3)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+3)) + ...
      EI(3)*Ux(3)*y(2*n+12+9*n+(3-1)*9+j)*cos(y(n+i)-y(n+3))+ ...
      EI(3)*Uy(3)*y(2*n+12+9*n+(3-1)*9+j)*sin(y(n+i)-y(n+3)) ); 
    end
end


for j=4:6  % alphas
    for i=1:3 % tube number
  if ~(i==j-3)
      dux(j,i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      EI(j-3)*Ux(j-3)*sin(y(n+i)-y(n+j-3))- EI(j-3)*Uy(j-3)*cos(y(n+i)-y(n+j-3)) );
      duy(j,i)=(1/(EI(1)+EI(2)+EI(3)))* (...
      EI(j-3)*Ux(j-3)*cos(y(n+i)-y(n+j-3))+ EI(j-3)*Uy(j-3)*sin(y(n+i)-y(n+j-3)) );
  end
    end
end

% with respect to uz0
% derivative of curvature gradient with respect to uz0
for j=7:9  
for i=1:3
dux(j,i)= (1/(EI(1)+EI(2)+EI(3)))* (...
     - EI(1)*Ux(1)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+1))+ ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+1))  + ...
       EI(1)*Ux(1)*y(2*n+12+9*n+(1-1)*9+j)*sin(y(n+i)-y(n+1))- ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(1-1)*9+j)*cos(y(n+i)-y(n+1))  + ...
       -EI(2)*Ux(2)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+2))+ ...
       EI(2)*Uy(2)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+2))  + ...
       EI(2)*Ux(2)*y(2*n+12+9*n+(2-1)*9+j)*sin(y(n+i)-y(n+2))- ...
       EI(2)*Uy(2)*y(2*n+12+9*n+(2-1)*9+j)*cos(y(n+i)-y(n+2))  + ...
       -EI(3)*Ux(3)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+3))+ ...
       EI(3)*Uy(3)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+3)) + ...
       EI(3)*Ux(3)*y(2*n+12+9*n+(3-1)*9+j)*sin(y(n+i)-y(n+3))- ...
       EI(3)*Uy(3)*y(2*n+12+9*n+(3-1)*9+j)*cos(y(n+i)-y(n+3)) );
    
duy(j,i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      -EI(1)*Ux(1)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+1))+ ...
      -EI(1)*Uy(1)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+1))  + ...
       EI(1)*Ux(1)*y(2*n+12+9*n+(1-1)*9+j)*cos(y(n+i)-y(n+1))+ ...
       EI(1)*Uy(1)*y(2*n+12+9*n+(1-1)*9+j)*sin(y(n+i)-y(n+1))  + ...
      -EI(2)*Ux(2)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+2))+ ...
      -EI(2)*Uy(2)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+2))  + ...
      EI(2)*Ux(2)*y(2*n+12+9*n+(2-1)*9+j)*cos(y(n+i)-y(n+2))+ ...
      EI(2)*Uy(2)*y(2*n+12+9*n+(2-1)*9+j)*sin(y(n+i)-y(n+2))  + ...
      -EI(3)*Ux(3)*y(2*n+12+9*n+(i-1)*9+j)*cos(y(n+i)-y(n+3))+ ...
      -EI(3)*Uy(3)*y(2*n+12+9*n+(i-1)*9+j)*sin(y(n+i)-y(n+3)) + ...
      EI(3)*Ux(3)*y(2*n+12+9*n+(3-1)*9+j)*cos(y(n+i)-y(n+3))+ ...
      EI(3)*Uy(3)*y(2*n+12+9*n+(3-1)*9+j)*sin(y(n+i)-y(n+3)) ); 
end
end


% odes for twist
for i=1:n      
    dydt(i)=  (  (EI(i))/(GJ(i))  ) * ( ux(i)* Uy(i) -  uy(i)* Ux(i) );  % ui_z
    dydt(n+i)=  y(i);   %alpha_i
    dydt(n+i)=  y(i);   %alpha_i
end

% derivative of z curvature  gradient with respect to inputs
for j=1:9
for i=1:n   
    if i==j
    duz(j,i)= ((dEI(i)*GJ(i)-dGJ(i)*EI(i))/(GJ(i)^2)) * (ux(i)* Uy(i)-uy(i)* Ux(i))+ ...
    (EI(i)/GJ(i)) * (ux(i)* dUy(i)-uy(i)* dUx(i)+dux(j,i)* Uy(i)-duy(j,i)* Ux(i));  
    else
    duz(j,i)= (EI(i)/GJ(i)) * ( dux(j,i)* Uy(i) -  duy(j,i)* Ux(i) ); 
    end
end
end
% 
% for j=4:9
%     for i=1:3
%         duz(j,i)= (EI(i)/GJ(i)) * ( dux(j,i)* Uy(i) -  duy(j,i)* Ux(i) );
%     end 
% end

dalpha= reshape(y(2*n+12+1: 2*n+12+9*n),[9,3]);

dydt(2*n+12+1: 2*n+12+9*n)= duz(:);
dydt(2*n+12+9*n+1: 2*n+12+9*n+9*n)= dalpha(:);

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

for i=1
    for j= 1:9
        
        du_hat=[0 -duz(j,i) duy(j,i) ; duz(j,i) 0 -dux(j,i) ; -duy(j,i) dux(j,i) 0 ];
        dR=reshape(y(2*n+12+9*n+9*n+(j-1)*9+1 : 2*n+12+9*n+9*n+(j-1)*9+9),[3,3])';
        ddRj= dR*u_hat+ R1* du_hat;
        dydt(2*n+12+9*n+9*n+9*9+(j-1)*3+1 : 2*n+12+9*n+9*n+9*9+(j-1)*3+3)=  ddRj*e3;  % gradient of dr
        ddRj=ddRj'; 
        dydt(2*n+12+9*n+9*n+(j-1)*9+1 : 2*n+12+9*n+9*n+(j-1)*9+9)=ddRj(:);
        
    end
end



dydt(2*n+1)=dr1(1);dydt(2*n+2)=dr1(2);dydt(2*n+3)=dr1(3);
dR=dR1'; 
dR=dR(:);
for i=4:12
   dydt(2*n+i)=dR(i-3);
end

end

%%
function [r1,r2,r3,Uz,angle,U_z,U_x,U_y] = moving_CTR_nmpc(q,uz_0,l,l_k)

param  % load tube parameters inside param.m file

% q1 o q3 are robot base movments, q3 to q6 are rbot base rotation angle.


B=q(1:3);  % length of tubes before template
%initial angles
alpha_1=q(4);
alpha=[q(4) q(5) q(6)];


% segmenting tubes  
% check all inputs must have n elements, n is number of tubes
[L,d_tip,EE,UUx,UUy] = segmenting(E,Ux,Uy,l,B,l_k);

SS=L;
for i=1:length(L)
    SS(i)=sum(L(1:i));
%     plot((B(1)+SS(i))*ones(1,10),1:10,'b' ,'LineWidth',2)
end

% S is segmented abssica of tube after template
 S=SS(SS+min(B)>0)+min(B);
 E=zeros(n,length(S)); Ux=E; Uy=E;
 for i=1:n
    E(i,:)=EE(i,SS+min(B)>0); Ux(i,:)=UUx(i,SS+min(B)>0); Uy(i,:)=UUy(i,SS+min(B)>0);
 end
 % each (i,j) element of above matrices correspond to the jth segment of
 % ith tube, 1st tube is the most inner

 %% Solving ode for shape

span=[0 S];       % vector of tube abssica starting at zero
Length=[]; r=[]; U_z=[]; U_x=[];U_y=[]; angle=[];  RR=[]; % solved length, curvatures, and twist angles
%U1_after=[0;0;0];             % 1st tube initial curvature at segment beginning
r0=[ 0 0 0]'; R0=[cos(alpha_1) sin(alpha_1) 0; -sin(alpha_1) cos(alpha_1) 0; 0 0 1];
R0=reshape(R0,[9,1]);
%alpha=alpha-B.*uz_0'; 

for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.0000001];
y0_1=[r0 ; R0];

y0_2=zeros(2*n,1);
y0_2(n+1:2*n)=alpha;
y0_2(1:n)=uz_0;

y_0=[y0_2; y0_1];

[s,y] = ode23(@(s,y) ode5(s,y,Ux(:,seg),Uy(:,seg),E(:,seg).*I',G.*J,n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i
shape=[y(:,2*n+1),y(:,2*n+2),y(:,2*n+3)];
Length=[Length; s];
r=[r; shape];
U_z=[U_z; y(:,1:n )];
r0=shape(end,:)';
R0=y(end,2*n+4:2*n+12)';
angle=[angle; y(:,1+n:2*n )];
RR=[RR; y(:,2*n+4:2*n+12 )];
uz_0=U_z(end,:)';
alpha=[y(end,n+1),y(end,n+2),y(end,n+3)]';


EI=E(:,seg).*I';GJ=G.*J;
Uxx=Ux(:,seg);Uyy=Uy(:,seg);
for k=1:length(s)  
ux(k)= (1/(EI(1)+EI(2)+EI(3)))* (...
      EI(1)*Uxx(1)  + ...
      EI(2)*Uxx(2)*cos(y(k,n+1)-y(k,n+2))+ EI(2)*Uyy(2)*sin(y(k,n+1)-y(k,n+2))  + ...
      EI(3)*Uxx(3)*cos(y(k,n+1)-y(k,n+3))+ EI(3)*Uyy(3)*sin(y(k,n+1)-y(k,n+3)) );
uy(k)= (1/(EI(1)+EI(2)+EI(3)))* (...
       EI(1)*Uyy(1) + ...
      -EI(2)*Uxx(2)*sin(y(k,n+1)-y(k,n+2))+ EI(2)*Uyy(2)*cos(y(k,n+1)-y(k,n+2))  + ...
      -EI(3)*Uxx(3)*sin(y(k,n+1)-y(k,n+3))+ EI(3)*Uyy(3)*cos(y(k,n+1)-y(k,n+3)) ); 
end

U_x=[U_x; ux'];
U_y=[U_y; uy'];

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
function dydt = ode5(~,y,Ux,Uy,EI,GJ,n)

dydt=zeros(2*n+12,1);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% second n elements of y are twist angles, alpha_i
% last 12 elements are r (position) and R (orientations), respectively


% calculating 1st tube's curvatures in x and y direction
ux=zeros(n,1);uy=zeros(n,1);

% calculating tube's curvatures in x and y direction
for i=1:n  
ux(i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      EI(1)*Ux(1)*cos(y(n+i)-y(n+1))+ EI(1)*Uy(1)*sin(y(n+i)-y(n+1))  + ...
      EI(2)*Ux(2)*cos(y(n+i)-y(n+2))+ EI(2)*Uy(2)*sin(y(n+i)-y(n+2))  + ...
      EI(3)*Ux(3)*cos(y(n+i)-y(n+3))+ EI(3)*Uy(3)*sin(y(n+i)-y(n+3)) );
uy(i)= (1/(EI(1)+EI(2)+EI(3)))* (...
      -EI(1)*Ux(1)*sin(y(n+i)-y(n+1))+ EI(1)*Uy(1)*cos(y(n+i)-y(n+1))  + ...
      -EI(2)*Ux(2)*sin(y(n+i)-y(n+2))+ EI(2)*Uy(2)*cos(y(n+i)-y(n+2))  + ...
      -EI(3)*Ux(3)*sin(y(n+i)-y(n+3))+ EI(3)*Uy(3)*cos(y(n+i)-y(n+3)) ); 
end

% odes for twist
for i=1:n      
    dydt(i)=  (  (EI(i))/(GJ(i))  ) * ( ux(i)* Uy(i) -  uy(i)* Ux(i) );  % ui_z
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

function [L,d1,E,Ux,Uy] = segmenting(E,Ux,Uy,l,B,l_k)

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
    I=E; Ux=E; Uy=E;
    return
end
end

EE=zeros(k,length(L));
 UUx=EE; UUy=EE;

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
UUx(i,b:c-1)=Ux(i);
UUy(i,b:c-1)=Uy(i);
end

l=L(~(L==0));  % get rid of zero lengthes
E=zeros(k,length(l)); Ux=E; Uy=E;
 for i=1:k
    E(i,:)=EE(i,~(L==0)); Ux(i,:)=UUx(i,~(L==0)); Uy(i,:)=UUy(i,~(L==0));
 end
L=L(~(L==0));

end

