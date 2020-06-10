% *************************************************************************
%             Shape Functions For Lagrange Elements (2D)
% *************************************************************************
% DESCRIPTION: Shape Functions for Lagrange Quad Elements such as QUAD4,
% QUAD9...
% MAT-files required: evaluateShapeFunction1D.m
%
% PARAMETERS:
%     p:    degree of the polynomial
%     x:    coordinates of a 2d point (x1, x2) in local coordinates system
%     f:    shape functions
%     df:   derivative of shape function
%
% EXAMPLES:
%     p = 2
%     x = [0.5, 0.5];
%     [shapeFunctions, derShapeFunctions] = evaluateShapeFunction2D(p, x);
%
% *************************************************************************
% Author: Hefeng Chen
% Email:  chhefeng@gmail.com
% December 2019; Last revision: 12-May-2020
% *************************************************************************

function [f, df] = evaluateShapeFunction2D(p, x)

% [NG, NDIM] = size(X);

NDIM = 2;
nBasisFunctions = p + 1;
nFunctions = nBasisFunctions ^ NDIM;

f  = zeros(1, nFunctions);
df = zeros(2, nFunctions);

[fx, dfx] = evaluateShapeFunction1D (p, x(1));
[fy, dfy] = evaluateShapeFunction1D (p, x(2));

t = 0;
for ny = 1: nBasisFunctions
    for nx = 1: nBasisFunctions
        
        t = t + 1;
        f (1, t) = fx(nx) * fy(ny);
        df(1, t) = dfx(:, nx) * fy(:, ny);
        df(2, t) = fx(:, nx)  * dfy(:, ny);
    end
end

end


