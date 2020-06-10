%% this script automates the grid convergence study for each method
% for plotting
nTicks = [1,2,3,4,5];

%% FTCS method
% initialise variables
ftcsT = zeros(1,5);
ftcsErr = zeros(1,5);
% calculate time taken and error for each grid resolution n=1,2,3,4,5
for n = 1:5
   [X,Y,U,t,Err] = ftcs(n);
   ftcsT(n) = t;
   ftcsErr(n) = Err;
   disp("FTCS simulation " + n + " has finished successfully");
end
disp("All FTCS simulations have now finished successfully");
% Plotting
% figure 1 is error magnitude relative to the analytical solution vs grid resolution
figure
plot(nTicks,ftcsErr,'b-s');
xticks([1 2 3 4 5]);
xticklabels({'n = 1','n = 2','n = 3','n = 4','n = 5'});
ylabel("Error Magnitude");
xlabel("Grid Resolution Level n");
title("Error relative to Analytical Solution vs Grid Resolution");
% figure 2 is error vs computation time
figure
plot(ftcsT,ftcsErr,'b-s');
ylabel("Error Magnitude");
xlabel("Computation Time (sec)");
title("Error relative to Analytical Solution vs Computation Time");

%% Central Differene method
% initialise variables
cdT = zeros(1,5);
cdErr = zeros(1,5);
% calculate time taken and error for each grid resolution n=1,2,3,4,5
for n = 1:5
   [X,Y,U,t,Err] = lpl_cd(n);
   cdT(n) = t;
   cdErr(n) = Err;
   disp("Iterative central difference simulation " + n + " has finished successfully");
end
disp("All Iterative central difference simulations have now finished successfully");
% Plotting
% figure 1 is error magnitude relative to the analytical solution vs grid resolution
figure
plot(nTicks,cdErr,'b-s');
xticks([1 2 3 4 5]);
xticklabels({'n = 1','n = 2','n = 3','n = 4','n = 5'});
ylabel("Error Magnitude");
xlabel("Grid Resolution Level n");
title("Error relative to Analytical Solution vs Grid Resolution");
% figure 2 is error vs computation time
figure
plot(cdT,cdErr,'b-s');
ylabel("Error Magnitude");
xlabel("Computation Time (sec)");
title("Error relative to Analytical Solution vs Computation Time");

%% Matrix method
% initialise variables
mtxT = zeros(1,5);
mtxErr = zeros(1,5);
% calculate time taken and error for each grid resolution n=1,2,3,4,5
for n = 1:5
   [X,Y,U,t,Err] = lpl_matrix(n);
   mtxT(n) = t;
   mtxErr(n) = Err;
   disp("Matrix central difference simulation " + n + " has finished successfully");
end
disp("All Matrix central difference simulations have now finished successfully");
% Plotting
% figure 1 is error magnitude relative to the analytical solution vs grid resolution
figure
plot(nTicks,mtxErr,'b-s');
xticks([1 2 3 4 5]);
xticklabels({'n = 1','n = 2','n = 3','n = 4','n = 5'});
ylabel("Error Magnitude");
xlabel("Grid Resolution Level n");
title("Error relative to Analytical Solution vs Grid Resolution");
% figure 2 is error vs computation time
figure
plot(mtxT,mtxErr,'b-s');
ylabel("Error Magnitude");
xlabel("Computation Time (sec)");
title("Error relative to Analytical Solution vs Computation Time");

% play chime when simulations finish
sound(sin(1:3000));
disp("All simulations have now finished successfully");