clear
close all
import casadi.*

%% Parameters:

% Space discretization
%dx = 0.2;
dx = 2/13;
nx = 2/dx-1;
sigma = 10.0;       
% Time discretization
nt = 200;          
% Final time
T = 0.1;                 

%% Building a grid
% Periodic boundary conditions in x, so nx+1 unknowns
xr = linspace(0, 2, nx+1);
% Dirichlet boundary conditions in y, so nx unknowns
yr_ = linspace(-1, 1, nx+2);
yr = yr_(2:nx+1);

%% Initial datum
Y0 = zeros(nx+1, nx);
H0 = zeros(nx+1, 1);
for i=1:nx+1
    H0(i, 1)= 2.5*xr(i)*(xr(i)-2);
    for j=2:nx
        %Y0(i,j) = 70*sin((xr(i)-1)*pi)*sin(yr(j)*pi);
        Y0(i,j) = 70*sin(xr(i)*pi)*sin(yr(j)*pi);
    end
end
% Vectorizing
W0 = zeros((nx+1)*nx+(nx+1), 1);
W0(1:(nx+1)*nx, 1) = Y0(:);
W0((nx+1)*nx+1:(nx+1)*nx+(nx+1),1) = H0;

%% Dynamics matrices A and B

% B. En gros ça veut dire que \omega 
% c'est une droite.
B = zeros((nx+1)*nx+(nx+1), (nx+1)*nx+(nx+1));
nu = 0.5*(nx+1)*nx;
B(0.25*(nx+1)*nx:0.75*(nx+1)*nx, 0.25*(nx+1)*nx:0.75*(nx+1)*nx) = eye(nu+1);


% We first define each block corresponding to the Laplacian with "pbc"
A0 = -4/dx^2*eye(nx+1)+1/dx^2*(diag(ones(nx, 1), 1)+diag(ones(nx, 1),-1));
A0(nx+1, 1) = 1/dx^2; 
A0(1, nx+1) = 1/dx^2;
% We also create the blocks of identities on the different "diagonals"
A1 = 1/dx^2*eye(nx+1);
% The other 3 blocks:
% Lowest on diagonal
A2 = -2*sigma/(2*dx^3)*eye(nx+1)+sigma/(2*dx^3)*(diag(ones(nx, 1), 1)+diag(ones(nx, 1),-1));
A2(nx+1, 1) = sigma/(2*dx^3);
A2(1, nx+1) = sigma/(2*dx^3);
% Residuals from equation for u:
A3 = -2*sigma/dx^4*eye(nx+1)+sigma/dx^4*(diag(ones(nx, 1), 1)+diag(ones(nx, 1),-1));
A3(nx+1, 1) = sigma/(dx^4);
A3(1, nx+1) = sigma/(dx^4);
% Free boundary residual of u: 
A4 = -1/(2*dx)*eye(nx+1);

% We now try to fill the matrix..
A = zeros((nx+1)*(nx+1), (nx+1)*(nx+1));
% First we fill the diagonal blocks with A0
for k = 1:nx
    l = (k-1)*(nx+1);
    A((l+1):(l+(nx+1)), (l+1):(l+(nx+1))) = A0;
end
% We now fill the off diagonal blocks
for k = 1:nx-1
    l = (k-1)*(nx+1);
    A((l+1):(l+(nx+1)), (l+(nx+1)+1):(l+2*(nx+1))) = A1;
end
 
for k = 1:nx-1
    l = (k-1)*(nx+1);
    A((l+(nx+1)+1):(l+2*(nx+1)), (l+1):(l+(nx+1))) = A1;
end
 
% We finally fill the 3 strange blocks
A(nx*(nx+1)+1:(nx+1)*(nx+1), nx*(nx+1)+1:(nx+1)*(nx+1)) = A2;
A((nx-1)*(nx+1)+1:nx*(nx+1), nx*(nx+1)+1:(nx+1)*(nx+1)) = A3;
A(nx*(nx+1)+1:(nx+1)*(nx+1), (nx-2)*(nx+1)+1:(nx-1)*(nx+1)) = A4;

%mu = 1/dx^2*0.001;
mu = 1;
A = sparse(mu*A);


%% Time-stepping: CFL condition to change 'nt'.
tline = linspace(0, T, nt+1); 
dt = tline(2)-tline(1);
CFL = mu*dt/dx^2;

%% ---- Input variables ---------
opti = casadi.Opti();               %% CasADi function
W = opti.variable((nx+1)^2, nt+1);      %% state trajectory
U = opti.variable((nx+1)^2, nt+1);        %% control

%% ---- Dynamic constraints --------
M = eye((nx+1)^2) - 0.5*dt*A;
L = eye((nx+1)^2) + 0.5*dt*A;
P = 0.5*dt*B;

for k=1:nt 
   % Crank-Nicolson method
   opti.subject_to(M*W(:,k+1) == L*W(:,k) + P*(U(:,k)+U(:,k+1)));
end

%% ---- State constraints --------
opti.subject_to(W(:,1)==W0);
opti.subject_to(W(:,nt+1)==0);

%% ---- Optimization objective  ----------
Cost = (dx*sum(sum(U.^2))*(T/nt));
opti.minimize(Cost); 

%% ---- initial guesses for solver ---
Y = zeros((nx+1)^2,nt+1);
opti.set_initial(W, Y);
opti.set_initial(U, 0);

%% ---- solve NLP ------
p_opts = struct('expand', true);
s_opts = struct('max_iter', 10000);     %% iteration limitation
opti.solver('ipopt', p_opts, s_opts);   %% set numerical backend

tic
sol = opti.solve();                     %% actual solve
Sol_u = sol.value(U);                   %% solved control function
Sol_x = sol.value(W);
time_axis = linspace(0, T, nt+1);

%% Plotting 
plot(time_axis, sum(Sol_u.^2)*dx, 'linewidth', 3, 'color', 'b')

ax = gca;
ax.LineWidth=1.5;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
ax.GridLineStyle = ':';
ax.MinorGridLineStyle = 'none';
exportgraphics(ax, 'stefan_hum_control.pdf', 'ContentType', 'vector')
toc

%% Plotting

% Splitting Y and H
Y = Sol_x(1:(nx+1)*nx, :);
H = Sol_x((nx+1)*nx+1:(nx+1)^2, :);

size(Y(:, 1))
% Plots
[Xr, Yr] = meshgrid(xr, yr_);

% higher res:
rez = 20;
xhr = linspace(0, 2, rez*(nx+1));
yhr_ = linspace(-1, 1, rez*(nx+2));
yhr = yhr_(2:rez*(nx+2)-1);
[Xhr, Yhr] = meshgrid(xhr, yhr_);

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'stefan-control.gif';
for t = 1:nt+1
    sol = reshape(Y(:,t), nx, nx+1);
    sol = [zeros(1, nx+1) ; sol ; transpose(sigma*A0*H(:, t))];
    solh = interp2(Xr, Yr, sol, Xhr, Yhr, 'cubic');
    ax = gca;
    z = surf(Xhr, Yhr, solh);
    %contourf(Xhr, Yhr, solh, 10);
    view(0, 90);
    colorbar;
    colormap("hot");
    lim = 1000;
    zlim([-lim lim])
    camlight
    lighting phong
%     z.AmbientStrength = 0.3;
%     z.DiffuseStrength = 0.8;
%     z.SpecularStrength = 0.9;
%     z.SpecularExponent = 25;
%     z.BackFaceLighting = 'unlit';
    shading interp
    caxis([-lim, lim])
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if t == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    if t==1 || t==10 || t==nt/4 || t==nt/2 || t==3*nt/4 || t==nt-10 || t==nt+1
        exportgraphics(ax, strcat('stefan_y_', num2str(t), '.pdf'), 'ContentType', 'vector');
    end
end

h1 = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename2 = 'stefan_h-control.gif';
for t = 1:nt+1
    hhr = interp1(xr, H(:, t), xhr, 'spline');
    plot(xhr, hhr, 'linewidth', 3, 'color', 'm');
    xlim([0 2]);
    ylim([-3.5 0.5]);
    ax = gca;
    ax.LineWidth=1.5;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    set(gca,'XMinorTick','on','YMinorTick','on')
    grid minor
    ax.GridLineStyle = ':';
    ax.MinorGridLineStyle = 'none';
    drawnow
    % Capture the plot as an image
    frame = getframe(h1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if t == 1
        imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename2,'gif','WriteMode','append');
    end
    if t==1 || t==10 || t==nt/4 || t==nt/2 || t==3*nt/4 || t==nt-10 || t==nt+1
    %if t==1 || t==floor((nt+1)/4) || t==floor(3*(nt+1)/4) || t==nt+1
        exportgraphics(ax, strcat('stefan_h_', num2str(t), '.pdf'), 'ContentType', 'vector');
    end
end