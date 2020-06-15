function [U,t,Err] = lpl_matrix(n)
%% Numerical approximation of Laplace equation
% n             - spatial step factor
% U             - solution values
% t, Err        - time and error for simulation
A = analyticalSol(n); % analyticalSol function for error calculations

% initialise stencil variables
L = 3;          % domain length
H = 1;          % domain height
ds = 0.1/2^n;   % spatial step size, depends on input n
T1 = -5;        % boundary temp 1
T2 = -10;       % boundary temp 2

% initialise solution variables
x = 0:ds:L;
y = 0:ds:H;
nx = length(x);
ny = length(y);
f1 = x;
f2 = x;
for pos = 2:nx-1
    f1(pos) = (-4*T1/L^2)*x(pos)*(x(pos)-L);
    f2(pos) = (-4*T2/L^2)*x(pos)*(x(pos)-L);
end

tic; % begin timing
% initialise output matrix
U = sparse(ny,nx);  % sparse matrix for larger n values
U(end,:) = f1;      % top boundary = f1
U(1,:) = f2;        % bottom boundary = f2
U(:,1) = 0;         % left boundary = 0
U(:,end) = 0;       % right boundary = 0

% initialise coefficients matrix M
NX = nx-2;              % solution points in x minus boundary conditions
NY = ny-2;              % solution points in y minus boundary conditions
size = NX*NY;           % unknown points
M = sparse(size,size);  % sparse matrix for larger n values
M(1:size+1:size^2) = -4;            % (i,j) diagonal
M(size+1:size+1:size^2) = 1;        % (i-1,j) diagonal
M(2:size+1:size^2-size) = 1;        % (i+1,j) diagonal
M(NX+1:size+1:size^2-NX*size) = 1;  % (i,j+1) diagonal
M(NX*size+1:size+1:size^2) = 1;     % (i,j-1) diagonal
M(size*(NX-1)+(NX+1):NX*size+NX:size^2) = 0; % i+1 at x = L boundary condition
M(size*NX+NX:NX*size+NX:size^2) = 0;         % i-1 at x = 0 boundary condition

% initialise solution vector Un
Unp1 = sparse(size,1); % solution vector

% initialise coefficients vector B
B = sparse(size,1); % boundary conditions are all zero except top and bottom
B(1:NX) = f2(2:nx-1); % bottom boundary
B(size-(NX-1):size) = f1(2:nx-1); % top boundary

% matrix multiply for solution vector
Unp1 = M\-B;

% supplant converged solution vector in output matrix
% note: reshape to get vector in matrix form, transpose to orientate correctly
U(2:end-1,2:end-1) = transpose(reshape(Unp1,[NX,NY]));
t = toc; % stop timing

%% error calculation - largest difference between two corresponding points in U and A
Err = max(abs(U(:)-A(:)));
end