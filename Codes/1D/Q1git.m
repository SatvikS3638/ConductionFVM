%Satvik Sarode, 2022A4PS0578P, 16 March
%Q1 - GIT & Analytical solution

clc
clear
close

% Material Properties

Ta = 100;
Tb = 500;
k = 1000;
A = 1e-2;
L = 0.5;
q = 1000;

%Analytical Solution
x_analytical = linspace(0, L, 100);
T_analytical = Ta + (Tb - Ta) * x_analytical / L + (q/(2*k)) * (x_analytical .* (L - x_analytical));
hold on;
figure(1);
plot(x_analytical, T_analytical, 'k--', 'LineWidth', 2, 'DisplayName', 'Analytical Solution');


%Grid Independence test variables
n_values = [5, 11, 21, 31, 41, 51, 61];
giterr = zeros(length(n_values), 1);

for m = 1:(length(n_values))
    % Discretisation of flow domain
    nn = n_values(m);
    dx = L/nn;          % Length of control volume domain
    M = zeros(nn,nn);
    S = zeros(nn, 1);
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
    
    for i = 1:nn
        M(i,i) = ap(i);
        S(i) = Su(i);
        if i ~= [1,nn]
            M(i,i+1) = -ae(i);
            M(i,i-1) = -aw(i);
        elseif i == 1
            M(i,i+1) = -ae(i);
        else 
            M(i,i-1) = -aw(i);
        end
    end
    
    T = linsolve(M,S);
    T = [Ta;T;Tb];

    %Grid Independence test
    midnod = (n_values(m)+1)/2;
    if m ~=1
        giterr(m) = (T(midnod, 1)-Tref(m-1))*100/Tref(m-1);
    end 
    Tref(m) = T(midnod, 1); 
    plot(X, T, '-o', 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(m)]);
end 

xlabel('x (m)'); ylabel('Temperature (°C)');
title('Steady-State 1D Heat Conduction with Heat Generation: CFD vs Analytical Solution');

figure(2);
plot(n_values(2:length(n_values)), giterr(2:length(n_values)));
xlabel("Mesh size");
ylabel("Relative error");
title("relative percentage error for variable mesh sizes");

