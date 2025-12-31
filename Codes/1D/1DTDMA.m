%Satvik Sarode, 2022A4PS0578P, 16 March
%Q1 - solution using TDMA

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
    Su(i) = q*A*dx;
    Sp(i) = 0;
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

b = sym(zeros(nn,1));
D = sym(zeros(nn,1));
a = sym(zeros(nn,1));
C = S;
A = sym(zeros(nn,1));
C1 = sym(zeros(nn,1));

for i = 1:nn
    if i == 1
        b(i) = 0;
        D(i) = M(i,i);
        a(i) = -M(i, i+1);
        A(i) = a(i)/(D(i));
        C1(i) = C(i)/D(i);
elseif i == nn
        b(i) = -M(i, i-1);
        D(i) = M(i,i);
        a(i) = 0;
        A(i) = a(i)/(D(i)-(b(i)*A(i-1)));
        C1(i) = (C(i)+(b(i)*C1(i-1)))/(D(i)-(b(i)*A(i-1)));
else
    b(i) = -M(i, i-1);
        D(i) = M(i,i);
        a(i) = -M(i, i+1);
        A(i) = a(i)/(D(i)-(b(i)*A(i-1)));
        C1(i) = (C(i)+(b(i)*C1(i-1)))/(D(i)-(b(i)*A(i-1)));
        end 
        end 
T = sym(zeros(nn, 1));    
for i=nn:-1:1
    if i == nn
        T(i) = C1(i);
    else
        T(i) = (A(i)*T(i+1))+C1(i);
    end
end
disp(T);
T = [Ta;T;Tb];

%Assuming q = 1000 W/m^3 as an example
q_value = 100000;
T_numeric = double(subs(T, q, q_value));

plot(X,T_numeric,"k+:");
xlabel("Distance(m)");
ylabel("Temperature(C°)");
title("Temperature Distribution Along the Thickness");
