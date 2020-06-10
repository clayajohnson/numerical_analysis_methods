function [X,Y,U,t,Err] = lpl_matrix(n)
%% Numerical approximation of 2D Laplace equation
% input:
% n             - spatial step factor
%
% output:
% X,Y,U         - meshgrid and solution values for plotting
% t             - time for error to converge
% L1,L2,LInf    - error norms for p = 1,2 and 10000

%% Analytical solution
% analyticalSol function is used to calculate the analytical solution for error norm calculations
A = analyticalSol(n);

%% Matrix operations Laplace central difference solution
% initialise stencil variables
L = 3;
H = 1;
ds = 0.1/2^n;   % spatial step size depends on input n
T1 = -5;        % boundary temp 1
T2 = -10;       % boundary temp 2

% define solution domain
x = 0:ds:L;
y = 0:ds:H;
nx = length(x);
ny = length(y);
f1 = x;
f2 = x;
[X,Y] = meshgrid(x,y);

% calculate boundary conditions
for pos = 2:nx-1
    f1(pos) = (-4*T1/L^2)*x(pos)*(x(pos)-L);
    f2(pos) = (-4*T2/L^2)*x(pos)*(x(pos)-L);
end

tic; % begin timing outside matrix initialisation
% initialise output matrix
U = sparse(ny,nx);
U(end,:) = f1;   % top boundary = f1
U(1,:) = f2;     % bottom boundary = f2
U(:,1) = 0;      % left boundary = 0
U(:,end) = 0;    % right boundary = 0

% initialise coefficients matrix M
NX = nx-2; % solution points in x minus boundary conditions
NY = ny-2; % solution points in y minus boundary conditions
size = NX*NY; % unknown points
M = sparse(size,size); % sparse matrix for larger n values
M(1:size+1:size^2) = -4; % (i,j) diagonal
M(size+1:size+1:size^2) = 1; % (i-1,j) diagonal
M(2:size+1:size^2-size) = 1; % (i+1,j) diagonal
M(NX+1:size+1:size^2-NX*size) = 1; % (i,j+1) diagonal
M(NX*size+1:size+1:size^2) = 1; % (i,j-1) diagonal
M(size*(NX-1)+(NX+1):NX*size+NX:size^2) = 0; % i+1 at x = L boundary condition
M(size*NX+NX:NX*size+NX:size^2) = 0; % i-1 at x = 0 boundary condition

% initialise solution vector Un
Unp1 = sparse(size,1); % solution vector

% initialise coefficients vector B
B = sparse(size,1); % boundary conditions are all zero except top and bottom
B(1:NX) = f2(2:nx-1); % bottom boundary
B(size-(NX-1):size) = f1(2:nx-1); % top boundary
B = -B;

% matrix multiply for new solution vector
Unp1 = M\B;

% supplant converged solution vector in output matrix
% note: reshape to get vector in matrix form, transpose to orientate correctly
U(2:end-1,2:end-1) = transpose(reshape(Unp1,[NX,NY]));

t = toc; % stop timing

%% error calculation
% the error is taken to be the largest difference between two corresponding points in U and A
Err = max(abs(U(:)-A(:)));
end