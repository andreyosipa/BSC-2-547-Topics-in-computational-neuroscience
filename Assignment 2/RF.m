function value = RF( x, y, sx, sy, phi, k )
% k in rad
value = (2*pi*sx*sy)^(-1) * exp(-x.^2/(2*sx^2) - y.^2/(2*sy^2)) * cos(k.*x./180.*pi - phi);
end

