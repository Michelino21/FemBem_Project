% Calculates the derivative of the kernel (1/r), which is part of the
% integrand for the coefficients of the H matrix.
%
%      met      Geometric element type (1 for linear elements, 2 for
%               quadratic).
%
%    (xi, yi)   The reference point.
%
%    (x1, y1) (x2, y2) (x3, y2)
%               The boundary points that identify the element. In case of
%               linear elements point (x2, y2) is the middle point between
%               (x1, y1) and (x3, y3).
%
%       t       The point in [-1, 1] where the function has to be
%               evaluated.

function [y] = intH(met, xi, yi, x1, y1, x2, y2, x3, y3, t)

    [rx, ry] = find_R(met, xi, yi, x1, y1, x2, y2, x3, y3, t);
    [nx, ny] = find_N(met, x1, y1, x2, y2, x3, y3, t);

    n = nx + 1i * ny;
    r = rx + 1i * ry;

    y = (1/(2*pi)) * cos(angle(n) - angle(r)) ./ abs(r);
    
    % theta = angle(n) - angle(r) is the angle between the radial direction
    % and the normal unit vector.

end