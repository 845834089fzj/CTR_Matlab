clear vars
clc
syms k1xy k1z k2xy k2z k3xy k3z tet1 a1 a2 a3 u1z u2z u3z u1x u2x u3x u1y u2y u3y U1x U2x U3x U1y U2y U3y U1z U2z U3z

% U is initial curvature of bent tube
U1=[U1x; U1y; U1z]; U2=[U2x; U2y; U2z]; U3=[U3x; U3y; U3z];
K1=diag([k1xy k1xy k1z]);
K2=diag([k2xy k2xy k2z]);
K3=diag([k3xy k3xy k3z]);
Rz1=[cos(a1) -sin(a1) 0; sin(a1) cos(a1) 0; 0 0 1];
Rz2=[cos(a2) -sin(a2) 0; sin(a2) cos(a2) 0; 0 0 1];
Rz3=[cos(a3) -sin(a3) 0; sin(a3) cos(a3) 0; 0 0 1];

da2=u2z-u1z;
da3=u3z-u1z;

uu1= (K1+K2+K3)^-1 * transpose(Rz1)*(Rz1*K1*U1 + Rz2*K2*U2 + Rz3*K3*U3);
uu2= (K1+K2+K3)^-1 * transpose(Rz2)*(Rz1*K1*U1 + Rz2*K2*U2 + Rz3*K3*U3);
uu3= (K1+K2+K3)^-1 * transpose(Rz3)*(Rz1*K1*U1 + Rz2*K2*U2 + Rz3*K3*U3);
u1x=simplify(uu1(1));
u1y=simplify(uu1(2));
u2x=simplify(uu2(1));
u2y=simplify(uu2(2));
u3x=simplify(uu3(1));
u3y=simplify(uu3(2));

duz1=(k1xy/k1z)*(u1x*U1y-u1y*U1x)
duz2=(k1xy/k1z)*(u1x*U1y-u1y*U1x)
duz3=(k1xy/k1z)*(u1x*U1y-u1y*U1x)


%% 

syms R11 R12 R13 R21 R22 R23 R31 R32 R33 u1 u2 u3 tetx(s) tety(s) tetz(s) dtetx dtety dtetz x

Rx=[1 0 0; 0 cos(tetx(s)) -sin(tetx(s)); 0 sin(tetx(s)) cos(tetx(s))];
Rz=[cos(tetz(s)) -sin(tetz(s)) 0; sin(tetz(s)) cos(tetz(s)) 0; 0 0 1];
Ry=[cos(tety(s)) 0 sin(tety(s)); 0 1 0; -sin(tety(s)) 0 cos(tety(s))];

R=Rx*Ry*Rz;


R=[R11 R12 R13; R21 R22 R23; R31 R32 R33];
e3=[0 0 1]';
%dr=R*e3;

u_hat=[ 0 -u3 u2; u3 0 -u1; -u2 u1 0];
%dR=R*u_hat;
%dR = diff(R,s);
%dR=subs(dR, [diff(tetx(s), s), diff(tety(s), s), diff(tetz(s), s)], [dtetx, dtety, dtetz]);

% 
% eqn = dR == R*u_hat;
% sol = solve(eqn,[dtetx, dtety, dtetz])

dr=expm(u_hat)*e3;

% u_hat_1=[ 0 0 0; u3 0 0; -u2 u1 0];
% u_hat_2=[  0 -u3 u2; 0 0 -u1; 0 0 0];
% 
% expm(u_hat_1)*expm(u_hat_2)

Rz*(eye(3,3)+(sin(x)/x)*u_hat+((1-cos(x))/x^2)*u_hat^2)*e3
