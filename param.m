% physical parameters of tubes

% Note: first tube is the most inner tube
n=3;                  % number of tubes
E1=30e9;             % Young's modulus
G1=11e9;             % Shear modulus
E2=30e9;             
G2=11e9;             
E3=30e9;
G3=11e9;
r_in1=0e-3;          % inner diameter tube 1
r_out1=1.25e-3;       % outer diameter tube 1
r_in2=1.75e-3;        % inner diameter tube 2 (outer)
r_out2=2.2e-3;      % outer diameter tube 2 (outer)
r_in3=2.5e-3;
r_out3=3e-5;

I1=(pi/4) * (r_out1^4-r_in1^4);   % 2nd moemnt of inertai tube 1
J1=(pi/2) * (r_out1^4-r_in1^4);   % polar moment of inertia tube 1
I2=(pi/4) * (r_out2^4-r_in2^4);   % 2nd moemnt of inertai tube 2
J2=(pi/2) * (r_out2^4-r_in2^4);   % polar moment of inertia tube 2
I3=(pi/4) * (r_out3^4-r_in3^4);   % 2nd moemnt of inertai tube 3
J3=(pi/2) * (r_out3^4-r_in3^4);   % polar moment of inertia tube 3

E=[E1 E2 E3]; I=[I1 I2 I3]; G=[G1 G2 G3]; J=[J1 J2 J3];
Ux=[-10 7 5];
Uy=[0 0 0];


%% Two tube

% n=2;                  % number of tubes
% E1=30e9;             % Young's modulus
% G1=11e9;             % Shear modulus
% E2=30e9;             
% G2=11e9;             
% E3=30e9;
% G3=11e9;
% r_in1=0e-3;          % inner diameter tube 1
% r_out1=1.25e-3;       % outer diameter tube 1
% r_in2=1.75e-3;        % inner diameter tube 2 (outer)
% r_out2=2.2e-3;      % outer diameter tube 2 (outer)
% r_in3=2.5e-3;
% r_out3=3e-5;
% 
% I1=(pi/4) * (r_out1^4-r_in1^4);   % 2nd moemnt of inertai tube 1
% J1=(pi/2) * (r_out1^4-r_in1^4);   % polar moment of inertia tube 1
% I2=(pi/4) * (r_out2^4-r_in2^4);   % 2nd moemnt of inertai tube 2
% J2=(pi/2) * (r_out2^4-r_in2^4);   % polar moment of inertia tube 2
% I3=(pi/4) * (r_out3^4-r_in3^4);   % 2nd moemnt of inertai tube 3
% J3=(pi/2) * (r_out3^4-r_in3^4);   % polar moment of inertia tube 3
% 
% E=[E1 E2]; I=[I1 I2]; G=[G1 G2]; J=[J1 J2];
% Ux=[10 5];
% Uy=[0 0];