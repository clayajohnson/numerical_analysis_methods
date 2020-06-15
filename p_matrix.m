function [X,Y,U,t] = p_matrix(n)
%% Numerical approximation of Poisson's equation
% n             - spatial step factor
% X,Y,U         - meshgrid and solution values for plotting
% t             - time for error to converge

% initialise stencil variables
L = 3;
H = 1;
ds = 1/(2*n^2);   % spatial step size depends on input n
P1 = 2;
P2 = -4;

% define solution domain
x = 0:ds:L;
y = 0:ds:H;
nx = length(x);
ny = length(y);
[X,Y] = meshgrid(x,y);

% initialise output matrix
U = sparse(ny,nx);
U(end,:) = 0;   % top boundary = 0
U(1,:) = 0;     % bottom boundary = 0
U(:,1) = 0;     % left boundary = 0
U(:,end) = 0;   % right boundary = 0

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

% initialise P vector
P = zeros(size,1);
p = 1;
for j = 1:NY
    for i = 1:NX
        P(p) = (P1*exp(-(((x(i)-1)^2)/0.1)-(((y(j)-0.5)^2)/0.05))+P2*exp(-(((x(i)-2)^2)/0.1)-(((y(j)-0.5)^2)/0.05)))*ds^2;
        p = p + 1;
    end
end
P = sparse(P); % convert to sparse

% initialise coefficients vector B
B = sparse(size,1); % boundary conditions are all zero
B = -B;

% matrix multiply for new solution vector
Unp1 = M\P + M\B;

% supplant converged solution vector in output matrix
% note: reshape to get vector in matrix form, transpose to orientate correctly
U(2:end-1,2:end-1) = transpose(reshape(Unp1,[NX,NY]));
end