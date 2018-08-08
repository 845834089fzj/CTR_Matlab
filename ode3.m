
function dydt = ode3(s,y,Uxx,Uyy,EI,GJ,tube_length,n)

dydt=zeros(2*n,1);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% second n elements of y are twist angles, alpha_i
% third n elements of y are curvatures along x
% fourth n elements of y are curvatures along y

% approximating physical parameters asfunction of tube length
ei=zeros(n,1); gj=ei; ux=ei; uy=ei;
for i=1:n
ei(i) = interp1(tube_length,EI(i,:),s); % Interpolate the data set (ft,f) at time t
gj(i) = interp1(tube_length,GJ(i,:),s); % Interpolate the data set (ft,f) at time t
ux(i) = interp1(tube_length,Uxx(i,:),s); % Interpolate the data set (ft,f) at time t
uy(i) = interp1(tube_length,Uyy(i,:),s); % Interpolate the data set (ft,f) at time t
end
EI=ei; GJ=gj; Ux=ux; Uy=uy;


% calculating summation of matrices
K=zeros(3,3);SUM=zeros(3,1);
for i=1:n
    k=diag([EI(i) EI(i) GJ(i)] );
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
    if GJ(i)==0
        GJ(i)=1;  % to avoid singularity when tube doesn't exist
    end
    dydt(i)=  (  (EI(i))/(GJ(i))  ) * ( ux(i)* Uy(i) -  uy(i)* Ux(i) );  % ui_z
    dydt(n+i)=  y(i)-y(1);   %alpha_i
end


end
