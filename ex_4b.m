clear all
close all
clc

syms theta dot_theta ddot_theta l x0 h m g real;

%Pc defitinion and its derivatives
pc=[l*cos(theta) + x0; 0; l*sin(theta) + h];
pc_d= [-l*sin(theta)*dot_theta; 0 ; l*cos(theta)*dot_theta];
pc_dd=[-l*cos(theta)*dot_theta^2-ddot_theta*l*sin(theta); 0; -l*sin(theta)*dot_theta^2+ddot_theta*l*cos(theta)];

%Skew Matrix Defitinion
S=[0 -1; 1  0];

%gravity vector definition
g0=[0; 0; -g];

% L matrix and its derivative computation
L=[0; m*l^2 * dot_theta; 0];
L_d=[0; m*l^2*ddot_theta; 0];

%Definition of vector pz
pz=pc(1:2)-(pc(3)/(pc_dd(3)-g0(3)))* (pc_dd(1:2)-g0(1:2))+(1/(m*(pc_dd(3)-g0(3))))*S*L_d(1:2);

pz=simplify(pz)