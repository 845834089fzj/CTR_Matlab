%% System of ODEs two tube

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
