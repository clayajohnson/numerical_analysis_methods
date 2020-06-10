uT = zeros(1,3);
for n = 1:3
    p = 1/(2*n^2);
    x = 0:p:3;
    y = 0:p:1;
    nx = length(x);
    ny = length(y);
    I = (ny*nx+1)/2;
    [X,Y,U] = p_matrix(n);
    figure
    [C,h] = contourf(X,Y,U);
    uT(n) = C(I);
    w = h.DataTipTemplate;
    w.DataTipRows(end) = dataTipTextRow('U',C);
    hold on;
    datatip(h,'DataIndex',I);
    xlabel("X");
    ylabel("Y    ",'Rotation', 0);
    title("Poisson's Equation for Grid Spacing " + p); 
    a = colorbar;
    ylabel(a,["Temperature " + char(176) + "C"], 'FontSize', 10, 'Rotation', 270);
    a.Label.Position(1) = 3;
    drawnow;
end

figure
plot(uT,'b-s');
xticks([1 2 3]);
xticklabels({1/(2^2),1/2^3,1/(2*3^2)});
xlabel("Grid Spacing");
ylabel("Temperature " + char(176) + "C");
title("Temperature at Centre of Domain for Varying Grid Spacing")
