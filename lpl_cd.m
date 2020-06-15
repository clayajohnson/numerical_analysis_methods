function [U,t,Err] = lpl_cd(n)
%% Numerical approximation of Laplace equation
% n             - spatial step factor
% U             - solution values
% t, Err        - time and error for simulation
A = analyticalSol(n); % analyticalSol function for error calculations

% initialise stencil variables
L = 3;          % domain length
H = 1;          % domain height
ds = 0.1/2^n;   % spatial step size depends on input n
T1 = -5;        % boundary temp 1
T2 = -10;       % boundary temp 2

% define solution domain
x = 0:ds:L;
y = 0:ds:H;
f1 = x;
f2 = x;
nx = length(x);
ny = length(y);

% calculate boundary conditions
for pos = 1:nx
    f1(pos) = (-4*T1/L^2)*x(pos)*(x(pos)-L);
    f2(pos) = (-4*T2/L^2)*x(pos)*(x(pos)-L);
end

% initialise solution matrix
Unp1 = zeros(ny,nx);
Unp1(end,:) = f1;   % top boundary y = H
Unp1(1,:) = f2;     % bottom boundary y = 0
Unp1(:,1) = 0;      % left boundary = 0
Unp1(:,end) = 0;    % right boundary = 0

% loop until error converges
err = 1;        % initial error value
tol = 1e-8;     % tolerance for error convergence
tic;            % begin timing
while  err > tol
    Un = Unp1; % update solution grid
    
    % loop over solution grid except for boundaries
    for j = 2:ny-1
        for i = 2:nx-1
            Unp1(j,i) = (1/4)*(Un(j+1,i) + Un(j-1,i) + Un(j,i+1) + Un(j,i-1));
        end
    end
    % calculte difference between iterations
    err = max(abs(Unp1(:)-Un(:)));
end
t = toc; % stop timing
U = Unp1; % Set output variable

%% error calculation - largest difference between two corresponding points in U and A
Err = max(abs(U(:)-A(:)));
end