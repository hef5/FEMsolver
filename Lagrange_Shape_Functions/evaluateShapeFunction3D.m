% *************************************************************************
%             Shape Functions For Lagrange Elements (3D)
% *************************************************************************
% DESCRIPTION: Shape Functions for Lagrange Cubic Elements such as HEX8,
% HEX27...
% MAT-files required: evaluateShapeFunction1D.m
%
% PARAMETERS:
%     p:    degree of the polynomial
%     x:    coordinates of a 3d point (x1) local coordinates system
%     f:    shape functions
%     df:   derivative of shape function
%
% EXAMPLES:
%     p = 2;
%     x = [0.2, 0.3, 0.3];
%     [shapeFunctions, derShapeFunctions] = evaluateShapeFunction1D(p, x);
%
% *************************************************************************
% Author: Hefeng Chen
% Email:  chhefeng@gmail.com
% December 2019; Last revision: 12-May-2020
% *************************************************************************

function [f, df] = evaluateShapeFunction3D(p, x)

NDIM = 3;
nBasisFunctions = p + 1;
nFunctions = nBasisFunctions ^ NDIM;

f  = zeros(1, nFunctions);
df = zeros(3, nFunctions);

[fx, dfx] = LagrangeBasis (x(1), p);
[fy, dfy] = LagrangeBasis (x(2), p);
[fz, dfz] = LagrangeBasis (x(3), p);

t = 0;
for nz = 1: nBasisFunctions
    for ny = 1: nBasisFunctions
        for nx = 1: nBasisFunctions
            t = t + 1;
            f (:, t) = fx(:, nx)  * fy(:, ny)  * fz(:, nz);
            df(1, t) = dfx(:, nx) * fy(:, ny)  * fz(:, nz);
            df(2, t) = fx(:, nx)  * dfy(:, ny) * fz(:, nz);
            df(3, t) = fx(:, nx)  * fy(:, ny)  * dfz(:, nz);
        end
    end
end

end
