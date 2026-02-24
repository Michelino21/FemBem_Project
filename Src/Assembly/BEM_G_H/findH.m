% Calculates the H matrix coefficient related to the given set of points.
%
%      met      Geometric element type (1 for linear elements, 2 for
%               quadratic).
%
%       c       Shape function type (2 is allowed only with quadratic
%               elements).
%
%    (xi, yi)   The reference point.
%
%    (x1, y1) (x2, y2) (x3, y2)
%               The boundary points that identify the element. In case of
%               linear elements point (x2, y2) is the middle point between
%               (x1, y1) and (x3, y3).

function [h_c] = findH(met, c, xi, yi, x1, y1, x2, y2, x3, y3)

    f = @(t) shape_fun(met, c, t)                                       ...
             .*                                                         ...
             intH(met, xi, yi, x1, y1, x2, y2, x3, y3, t)               ...
             .*                                                         ...
             abs_J(x1, y1, x2, y2, x3, y3, t);

    h_c = quadgk(f, -1, +1);

end