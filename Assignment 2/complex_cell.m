function [ out ] = complex_cell( image, range )
phi1 = pi/2;
phi2 = 0;

adj_1 = 50/4.6790e-07;
adj_2 = 50/0.9999;

nonlinearity = @(x) x^2;

out1 = simple_cell(image, range, phi1, nonlinearity, adj_1);
out2 = simple_cell(image, range, phi2, nonlinearity, adj_2);

out = out1 + out2;

end