% Calculates the kernel (ln r), which is part of the integrand for the
% the coefficients of the G matrix.
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

function [y] = intG(met, xi, yi, x1, y1, x2, y2, x3, y3, t)

    [rx, ry] = find_R(met, xi, yi, x1, y1, x2, y2, x3, y3, t);

    r = rx + 1i * ry;
    y = (1/(2*pi)) * log(abs(r));

end