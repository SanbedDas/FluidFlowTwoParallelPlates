clc
clear
%% Parameters
Lx = 1;      % Plate Length in x direction
Ly = 1;      % Gap between plates in y direction
nx = 20;     %  Grid points in x direction
ny = 20;     % Grid points in y direction
dx = Lx/(nx-1); % Grid spacing in x direction
dy = Ly/(ny-1); % Grid spacing in y direction
nt = 100; % Number of time steps
nu = 0.5; % Kinematic viscosity 
dt=0.001;
%% Initial Conditions

u = zeros(ny,nx); % x-velocity
v = zeros(ny,nx); % y-velocity
p = 10*ones(ny,nx); % pressure

%% Boundary Conditions
% Applying top wall boundary condition
u(1,:) = 10; 
for i=1:2:length(u)
    p(i,:)=15;
end

%% SIMPLE Algorithm
% Solve using semi-implicit pressure linked equations method
for n = 1:nt
    % Calculate intermediate velocity field
    u_star = u;
    v_star = v;
    for i = 2:nx-1
        for j = 2:ny-1
            u_star(j,i) = u(j,i) + dt*(-(u(j,i)*(u(j,i)-u(j,i-1))/dx)-(v(j,i)*(u(j,i)-u(j-1,i))/dy)...
                                       +(nu/dx^2)*(u(j,i+1)-2*u(j,i)+u(j,i-1))...
                                       +(nu/dy^2)*(u(j+1,i)-2*u(j,i)+u(j-1,i)));
            v_star(j,i) = v(j,i) + dt*(-(u(j,i)*(v(j,i)-v(j,i-1))/dx)-(v(j,i)*(v(j,i)-v(j-1,i))/dy)...
                                       +(nu/dx^2)*(v(j,i+1)-2*v(j,i)+v(j,i-1))...
                                       +(nu/dy^2)*(v(j+1,i)-2*v(j,i)+v(j-1,i)));
            u(:,nx)=u(:,nx-1);  
        end
    end
    
    rhs = zeros(ny,nx);
    for i = 2:nx-1
        for j = 2:ny-1
            rhs(j,i) = (1/dt)*((u_star(j,i)-u_star(j,i-1))/dx+(v_star(j,i)-v_star(j-1,i))/dy);
        end
    end
    p_new =p;
%     zeros(ny,nx);
    max_err = 1;
    while max_err > 1e-3
    p_old = p_new;
    for i = 2:nx-1
        for j = 2:ny-1
            p_new(j,i) = ((p_old(j,i+1)+p_old(j,i-1))*dy^2+(p_old(j+1,i)+p_old(j-1,i))*dx^2-rhs(j,i)*dx^2*dy^2)/(2*(dx^2+dy^2));
        end
    end
    max_err = max(abs(p_new(:)-p_old(:)));
    end  
    % Correct velocity field with pressure
    for i = 2:nx-1
        for j = 2:ny-1
            u(j,i) = u_star(j,i)-(dt/dx)*(p_new(j,i+1)-p_new(j,i));
            v(j,i) = v_star(j,i)-(dt/dy)*(p_new(j+1,i)-p_new(j,i));
        end
    end    
end
%% Plotting the velocity and pressure fields
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);
figure();
contourf(X,Y,flipud(u),20,'LineColor','none');
xlabel('x');
ylabel('y');
title(sprintf('Velocity field at t=%.2f',nt));
colorbar;

figure();
pressure=p_new+p_old;
contourf(X,Y,flipud(pressure),20,'LineColor','none');
xlabel('x');
ylabel('y');
title(sprintf('Pressure field at t=%.2f',nt));
colorbar;
% figure();