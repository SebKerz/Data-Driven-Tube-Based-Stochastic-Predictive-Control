function [A, B, Bd, T] = doubleMassOscillator()
% This is used to retrieve the discrete time model of the double
% mass oscillator, as well as an LQR feedback gain.
%% Parameters
J1 = 5e-2; %kg m^2
J2 = 5e-2; %kg m^2
d1 = 1e-2; %Nm/(rad/s)
d2 = 1e-2; %Nm/(rad/s)
c1 = 1; %N/rad

%New to test, 26.02.2023
J1 = 10e-2; %kg m^2
J2 = 10e-2; %kg m^2
d1 = 1e-1; %Nm/(rad/s)
d2 = 1e-1; %Nm/(rad/s)
c1 = 1; %N/rad


%% System matrices (continuous time)
A_con = [0, 0, 1, 0;
         0, 0, 0, 1;
         -c1/J1, c1/J1, -d1/J1, 0;
         c1/J2, -c1/J2, 0, -d2/J2];

B_con = [0; 
         0;  
         1/J1; 
         0];

Bd_con = [0; 
          0;  
          0; 
          1/J2];

%% Discretization
%Sampling time:
T = 0.1;

A = expm(A_con.*T);

B = integral(@(t) expm(A_con.*t)*B_con, 0, T, 'ArrayValued', 1);

Bd = integral(@(t) expm(A_con.*t)*Bd_con, 0, T, 'ArrayValued', 1);

