function [U] = analyticalSol(n)
%% Analytical solution to a 2D heat problem
% n             - spatial step factor
% U             - analytical solution for 100 non-zero fourier terms

% solution parameters
ds = 0.1/2^n;   % spatial grid size
L = 3;          % x domain is (0,L)
H = 1;          % y domain is (0,H)
T1 = -5;        % temperature for u1
T2 = -10;       % temperature for u2

% initialise domain
x = 0:ds:L;  
y = 0:ds:H;
nx = length(x);
ny = length(y);

for pos = 1:nx
    f1(pos) = (-4*T1/L^2)*x(pos)*(x(pos)-L);
    f2(pos) = (-4*T2/L^2)*x(pos)*(x(pos)-L);
end

% initialise solution output
U = zeros(ny,nx);
U(end,:) = f1;
U(1,:) = f2;

%time the computation
tic;
% run the calculation for 100 non-zero fourier terms
for nfs = 1:2:200
    bN = (2*cos(nfs*pi) - 2)/((nfs^3)*(pi^3)); % calculate fourier coefficient
    
    % loop over each x point
    for i = 2:nx-1
        xTerm = sin(nfs*pi*x(i)/L); % calculate current x component

        % loop over each y point
        for j = 2:ny-1
            u1 = (-8)*T1*bN*xTerm*sinh(nfs*pi*y(j)/L)/sinh(nfs*pi*H/L); % calculate one unique point with u1
            u2 = (-8)*T2*bN*xTerm*sinh(nfs*pi*(y(j)-H)/L)/sinh(((-nfs)*pi*H)/L); % calculate one unique point with u2
            U(j,i) = U(j,i) + u1 + u2; % add solutions from u1 and u2
        end
    end
end
% finish timing
t = toc;
end