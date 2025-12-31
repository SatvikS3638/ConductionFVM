%Satvik Sarode, 2022A4PS0578P, 16 March
%Modified program to study heat generation term
clc
clear
close

% Material Properties

Ta = 100;
Tb = 500;
k = 1000;
A = 1e-2;
L = 0.5;

syms q;
% Discretisation of flow domain

nn = 5;             % Total number of Nodes
dx = L/nn;          % Length of control volume domain
X = [0,dx/2:dx:L-dx/2,L];

% Solving Governing Equations for Non Extermal Nodes

for i = 2:nn-1
    aw(i) = k*A/dx;
    ae(i) = k*A/dx;
    Sp(i) = 0;
    Su(i) = q*A*dx;
    ap(i) = ae(i) + aw(i) - Sp(i);
end

for i = 1
    aw(i) = 0;
    ae(i) = k*A/dx;
    Sp(i) = -2*k*A/dx;
    Su(i) = 2*k*A*Ta/dx + q*A*dx;
    ap(i) = ae(i) + aw(i) - Sp(i);
end

for i = nn
    ae(i) = 0;
    aw(i) = k*A/dx;
    Su(i) = 2*k*A*Tb/dx +q*A*dx;
    Sp(i) = -2*k*A/dx;
    ap(i) = ae(i) + aw(i) - Sp(i);
end

% Coefficient Matrix
M = sym(zeros(nn,nn));
S = sym(zeros(nn, 1));

for i = 1:nn
    M(i,i) = ap(i);
    S(i) = Su(i);
    if i ~= [1,5]
        M(i,i+1) = -ae(i);
        M(i,i-1) = -aw(i);
    elseif i == 1
        M(i,i+1) = -ae(i);
    else 
        M(i,i-1) = -aw(i);
    end
end

fprintf("\nThe coefficient matrxix is \n");
disp(M);

T = linsolve(M,S);
fprintf("\nThe Temperature Distribution Vector is \n");
disp(T);
T = [Ta;T;Tb];
 
%Assuming q = 1000 W/m^3 as an example
q_value = 1000;
T_numeric = double(subs(T, q, q_value));
disp(T_numeric);
plot(X,T_numeric,"k+:");
xlabel("Distance(m)");
ylabel("Temperature(C°)");
title("Temperature Distribution Along the Thickness");

