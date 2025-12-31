%Satvik Sarode, 2022A4PS0578P, 16 March
%Q2 - Grid indeendece test for two heat sources


clc
clear
close

% Temperatures at the boundaries

Tl = 100;
Tr = 100;
Tt = 200;
Tb = 100;
q = 1000; %Input q value 

% Material Properties

k = 1;
t = 1;
L = 1;
B = 1;

% Grid Indeendence Test variables
n_values = [5, 7, 9, 11, 21, 31, 41, 51, 61, 71, 81, 91];
giterr = zeros(length(n_values), 1);

for m = 1:(length(n_values))
    % Discretisation of flow domain

nx = n_values(m);                           % Control volumes in x direction
ny = n_values(m);                           % Control volumes in y direction
nnx = nx + 2;                     % Total number of Nodes in x direction
nny = ny + 2 ;                    % Total number of Nodes in y direction
dx = L/nx;                        % Length of control volume in x direction
dy = B/ny;                        % Length of control volume in y direction
Ae = dy*t;
Aw = dy*t;
An = dx*t;
As = dx*t;
S = zeros(nx,ny);                 % Source Term Vector
X = [0,dx/2:dx:L-dx/2,L];         % X coordinates
Y = [0,dy/2:dy:L-dy/2,B];         % Y coordinates
phi = zeros(nx+2,ny+2);
phi_old = zeros(nx+2,ny+2);
err = zeros(nx+2,ny+2);
rms = 1;

% Solving Governing Equations for Non Extermal Nodes
% IMPORTANT!! - i corresponds to row transversal and j corresponds to column transversal%

% Interior Nodes
for i = 2:ny-1
    for j = 2:nx-1
        if (i == 3 && j == 4) || (i == 5 && j ==4)
            aw(i,j) = k*Aw/dx;
            ae(i,j) = k*Ae/dx;
            an(i,j) = k*An/dy;
            as(i,j) = k*As/dy;
            Sp(i,j) = 0;
            Su(i,j) = q*0.25;
            ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
        elseif (i == 2 && j == 4) || (i == 6 && j ==4)
            aw(i,j) = k*Aw/dx;
            ae(i,j) = k*Ae/dx;
            an(i,j) = k*An/dy;
            as(i,j) = k*As/dy;
            Sp(i,j) = 0;
            Su(i,j) = q*0.75;
            ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
        else 
        aw(i,j) = k*Aw/dx;
        ae(i,j) = k*Ae/dx;
        an(i,j) = k*An/dy;
        as(i,j) = k*As/dy;
        Sp(i,j) = 0;
        Su(i,j) = 0;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
        end 
    end
end
% Bottom Boundary Nodes
for i = ny
    for j = 2:nx-1
        aw(i,j) = k*Aw/dx;
        ae(i,j) = k*Ae/dx;
        an(i,j) = k*An/dy;
        as(i,j) = 0;
        Sp(i,j) = -2*k*As/dy;
        Su(i,j) = 2*k*As*Tb/dy;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Top Boundary Nodes
for i = 1
    for j = 2:nx-1
        aw(i,j) = k*Aw/dx;
        ae(i,j) = k*Ae/dx;
        an(i,j) = 0;
        as(i,j) = k*As/dy;
        Sp(i,j) = -2*k*An/dy;
        Su(i,j) = 2*k*An*Tt/dy;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Left Boundary Nodes
for j = 1
    for i = 2:ny-1
        aw(i,j) = 0;
        ae(i,j) = k*Ae/dx;
        an(i,j) = k*An/dy;
        as(i,j) = k*As/dy;
        Sp(i,j) = -2*k*Aw/dx;
        Su(i,j) = 2*k*Aw*Tl/dx;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Right Boundary Nodes
for j = nx
    for i = 2:ny-1
        aw(i,j) = k*Aw/dx;
        ae(i,j) = 0;
        an(i,j) = k*An/dy;
        as(i,j) = k*As/dy;
        Sp(i,j) = -2*k*Ae/dx;
        Su(i,j) = 2*k*Ae*Tr/dx;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Top Left node
for i = 1
    for j = 1
        aw(i,j) = 0;
        ae(i,j) = k*Ae/dx;
        an(i,j) = 0;
        as(i,j) = k*As/dy;
        Sp(i,j) = -2*k*Aw/dx - 2*k*An/dy;
        Su(i,j) = 2*k*Aw*Tl/dx + 2*k*An*Tt/dy;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Bottom Left node
for i = ny
    for j = 1
        aw(i,j) = 0;
        ae(i,j) = k*Ae/dx;
        an(i,j) = k*An/dy;
        as(i,j) = 0;
        Sp(i,j) = -2*k*Aw/dx - 2*k*As/dy;
        Su(i,j) = 2*k*Aw*Tl/dx + 2*k*As*Tb/dy;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Bottom Right node
for i = nx
    for j = ny
        aw(i,j) = k*Aw/dx;
        ae(i,j) = 0;
        an(i,j) = k*An/dy;
        as(i,j) = 0;
        Sp(i,j) = -2*k*Ae/dx - 2*k*As/dy;
        Su(i,j) = 2*k*Ae*Tr/dx + 2*k*As*Tb/dy;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Top Right node
for i = 1
    for j = ny
        aw(i,j) = k*Aw/dx;
        ae(i,j) = 0;
        an(i,j) = 0;
        as(i,j) = k*As/dy;
        Sp(i,j) = -2*k*Ae/dx - 2*k*An/dy;
        Su(i,j) = 2*k*Ae*Tr/dx + 2*k*An*Tt/dy;
        ap(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) - Sp(i,j);
    end
end

% Boundary Conditions
phi(1,2:nx+1) = Tt;
phi(ny+2,2:nx+1) = Tb;
phi(2:ny+1,1) = Tl;
phi(2:ny+1,nx+2) = Tr;

phi(1,1) = (Tt+Tl)/2;
phi(1,nx+2) = (Tt+Tr)/2;
phi(ny+2,1) = (Tb+Tl)/2;
phi(ny+2,nx+2) = (Tb+Tr)/2;

iter = 0;
errmax = 1e-4;   % Convergence criteria
rms = 1;
a = 0.7;         % Under relaxed case                       

while(iter<20000 && rms>errmax)
    iter = iter+1;
    phi_old = phi;
    for i = 2:ny+1
        for j = 2:nx+1
            phi(i,j) = (1-a)*phi(i,j) + (a/ap(i-1,j-1))*(aw(i-1,j-1)*phi(i-1,j) + ae(i-1,j-1)*phi(i+1,j) + an(i-1,j-1)*phi(i,j+1) + as(i-1,j-1)*phi(i,j-1) + Su(i-1,j-1));
        end
    end
    err = abs(phi - phi_old);
    rms = sqrt(sumsqr(err));
end
fprintf("\nThe Solution Converged at - %d iterations \n",iter);
fprintf("\nThe Temperature Distribution of the 2D Plate is - \n");
%disp(phi);
midnod = (n_values(m)+3)/2;

if m ~=1
    giterr(m) = (phi(midnod, midnod)-Tref(m-1))*100/Tref(m-1);
end 

Tref(m) = phi(midnod, midnod); 
if n_values(m) == 5
    % Post Processing the results
    X = fliplr(X);
    Y = fliplr(Y);
    [Xa,Ya] = meshgrid(X,Y);
    figure(1);
    [C,h]= contourf(Xa,Ya,phi);
    clabel(C,h);
    title('Contour plot of the Isotherms'),xlabel('X Coordinates'),ylabel('Y Coordinates'),colorbar;
    pos1 = get(gcf,'Position');                     % get position of Figure(1) 
    set(gcf,'Position', pos1 - [pos1(3)/2,0,0,0])   % Shift position of Figure(1)
    
    figure(2);
    mesh(Xa,Ya,phi);title('Mesh surface plot'),xlabel('x'),ylabel('y'),colorbar;
    pos2 = get(gcf,'Position');                     % get position of Figure(2) 
    set(gcf,'Position', pos2 + [pos1(3)/2,0,0,0])   % Shift position of Figure(2)
end
figure(3);
plot(n_values(2:length(n_values)), giterr(2:length(n_values)));
xlabel("Mesh size");
ylabel("Relative error");
end 

