nTicks = [1,2,3,4,5];
times = zeros(1,5);

for n = 1:5
    times(n) = analyticalSol(n);
    disp("Simulation for n = " + n + " has finished successfully");
end

figure
plot(nTicks,times,'b-s');
xticks([1 2 3 4 5]);
xticklabels({'n = 1','n = 2','n = 3','n = 4','n = 5'});
xlabel("Grid Resolution Level n");
ylabel("Computation Time (sec)");
title("Computation Time vs Grid Resolution");

disp("All done now thank fuck");
sound(sin(1:3000));