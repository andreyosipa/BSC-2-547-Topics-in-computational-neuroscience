function [ out ] = simple_cell( image, range, phi, nonlinearity, const )
%Assuming range is values of X(Y has same)
%Nonlinearity is some function, e.g. x^2

sx = 1;
sy = 2;
k = 1/0.56/180*pi; %convert k to rad

[X,Y] = meshgrid(range);
rf_values = zeros(length(range));
%evaluate RF:
for ix=1:length(range)
   for iy=1:length(range)
      rf_values(ix,iy) = RF( X(ix, iy), Y(ix, iy), sx, sy, phi, k);
   end
end
%adj constant, based on the values of RF.
adj = sum(sum(abs(rf_values)));

rf_values = rf_values ./ adj;

out = const * nonlinearity(sum(sum(rf_values .* image)));

end

