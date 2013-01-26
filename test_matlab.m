
grid_min = -2;
grid_max = 2;
grid_size = 10;
xi_min = grid_min * 1.2;
xi_max = grid_max * 1.2
xi_size = grid_size * 10;

% 1D
grid_1 = linspace(grid_min, grid_max, grid_size);
xi_grid_1 = linspace(xi_min, xi_max, xi_size);
X1 = grid_1;
f = sin(X1);
[xi_1] = ndgrid(xi_grid_1);
yi_1 = interp1(X1,f,xi_1); 
figure();
plot(X1,f,'o',xi_1,yi_1,'x')

% yi2 = linterp_matlab(x, y, xi);
% figure();
% plot(x,y,'o',xi,yi2)