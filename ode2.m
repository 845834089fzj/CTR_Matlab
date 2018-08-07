%% System of ODEs two tube

function dydt = ode2(s,y,U_x,U_y,U_z,Length)

e3=[0 0 1]';              


% approximating curvature of 1st tube as a function of tube length

ux = interp1(Length,U_x(:,1),s); 
uy = interp1(Length,U_y(:,1),s); 
uz = interp1(Length,U_z(:,1),s); 



% y(1) to y(3) are position of point materials
%r1=[y(1); y(2); y(3)];
% y(4) to y(12) are rotation matrix elements
R1=[y(4) y(5) y(6);y(7) y(8) y(9);y(10) y(11) y(12)];


u_hat=[0 -uz uy ; uz 0 -ux ; -uy ux 0 ];


% odes
dr1 = R1*e3;
dR1=R1*u_hat;


dydt(1)=dr1(1);dydt(2)=dr1(2);dydt(3)=dr1(3);
dR=dR1'; 
dR=dR(:);
for i=4:12
   dydt(i)=dR(i-3);
end


dydt=dydt';end
