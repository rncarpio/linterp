% change this to point to boost.multi_array
% mex -IC:/boost/boost_1_49_0 linterp_matlab.cpp

grid_min = -2;
grid_max = 2;
grid_size = 10;
xi_min = grid_min * 1.2;
xi_max = grid_max * 1.2;
xi_size = grid_size * 10;

% function to interpolate
sin_sum_1 = @(x) sin(x);
sin_sum_2 = @(x1, x2) sin(x1+x2);

% max error, error sum of squares
err_ss = @(true_y, y) [max(abs(y-true_y)), sum((y-true_y) .* (y-true_y))];

% 1D
grid_1 = linspace(grid_min, grid_max, grid_size);		% original grid
xi_grid_1 = linspace(xi_min, xi_max, xi_size);			% grid of points to interpolate on
f = sin_sum_1(grid_1);									% f evaluated on original grid
xi_mesh_1 = xi_grid_1;
true_y = sin_sum_1(xi_mesh_1);

% matlab's interp1
yi_1 = interp1(grid_1,f,xi_mesh_1); 					% interpolated value
not_nans = not(isnan(yi_1));
disp('interp1: max err, err SS');
disp(err_ss(true_y(not_nans), yi_1(not_nans)));
figure('Name', 'interp1');
plot(grid_1,f,'o', xi_grid_1,yi_1,'x')

% linterp
yi_2 = linterp_matlab(grid_1,f,xi_grid_1); 				% interpolated value
disp('linterp: max err, err SS');
disp(err_ss(true_y(not_nans), yi_2(not_nans)));
figure('Name', 'linterp 1d');
plot(grid_1,f,'o', xi_grid_1,yi_2,'x')

% 2D
grid_1 = linspace(grid_min, grid_max, grid_size);		% original grid
grid_2 = linspace(grid_min, grid_max, grid_size);
xi_grid_1 = linspace(xi_min, xi_max, xi_size);			% grid of points to interpolate on
xi_grid_2 = linspace(xi_min, xi_max, xi_size);
[grid_mesh_1, grid_mesh_2] = ndgrid(grid_1, grid_2);
[xi_mesh_1, xi_mesh_2] = ndgrid(xi_grid_1, xi_grid_2);
f = sin_sum_2(grid_mesh_1, grid_mesh_2);				% f evaluated on original grid
true_y = sin_sum_2(xi_mesh_1, xi_mesh_2);
figure('Name', 'true function');
mesh(grid_1, grid_2, f);

% matlab's griddata 2d
disp('griddata:');
tic;
yi_1 = griddata(grid_mesh_1, grid_mesh_2, f, xi_mesh_1, xi_mesh_2); 				% interpolated value
toc;
not_nans = not(isnan(yi_1));
disp(err_ss(true_y(not_nans), yi_1(not_nans)));
figure('Name', 'griddata 2d');
mesh(xi_grid_1, xi_grid_2, yi_1);
err = yi_1 - true_y;
figure('Name', 'griddata 2d errors');
mesh(xi_grid_1, xi_grid_2, err);

% linterp 2d
disp('linterp 2d:');
tic;
yi_2 = reshape(linterp_matlab(grid_1, grid_2, f, xi_mesh_1, xi_mesh_2), size(true_y)); 							% interpolated value
toc;
disp(err_ss(true_y(not_nans), yi_2(not_nans)));
figure('Name', 'linterp 2d');
mesh(xi_grid_1, xi_grid_2, yi_2);
err = yi_2 - true_y;
err(isnan(yi_1)) = NaN;
figure('Name', 'linterp 2d errors');
mesh(xi_grid_1, xi_grid_2, err);

