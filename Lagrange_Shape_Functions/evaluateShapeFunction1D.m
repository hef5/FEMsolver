% *************************************************************************
%             Shape Functions For Lagrange Elements (1D)
% *************************************************************************
% DESCRIPTION: Shape Functions for Lagrange Quad Elements such as BAR2,
% BAR3...
% MAT-files required: none
% 
% PARAMETERS:
%     p:    degree of the polynomial       
%     x:    coordinates of a 1d point (x1) in local coordinates system                   
%     f:    shape functions 
%     df:   derivative of shape function
%
% EXAMPLES:
%     p = 2;
%     x = 0.5;
%     [shapeFunctions, derShapeFunctions] = evaluateShapeFunction1D(p, x);
%
% *************************************************************************
% Author: Hefeng Chen
% Email:  chhefeng@gmail.com
% December 2019; Last revision: 12-May-2020
% *************************************************************************

function [shapeFunctions, derShapeFunctions] = evaluateShapeFunction1D(p, x)

    % check input parameters
    %     if (p < 1)
    %         error('Polynomial degree p is a integer and greater than 1');
    %     end
    %     if (x > 1 || x < -1)
    %         error('The range of coordinates value muss be [-1, 1]');
    %     end
    
    % number of shape function in 1D
    nShapeFunctions = p + 1;
    
    % the natural coordinate of node
    delta = 2 / p;
    
    xi = -1: delta: 1;
    
    % shape functions
    shapeFunctions = ones(1, nShapeFunctions);
    
    for i = 1:nShapeFunctions
    
        for k = 1:nShapeFunctions
            if (k ~= i)
              shapeFunctions(i) = shapeFunctions(i) * (xi(k)-x)/(xi(k)-xi(i));
            end
        end
        
    end
    
    % derivative of shape functions
    derShapeFunctions(nShapeFunctions) = 0;
    for i = 1:nShapeFunctions
             
        for l = 1:nShapeFunctions
            if (l ~= i)
                
                coefficient = 1;
                for k = 1:nShapeFunctions                   
                    if (k ~= i) && (k ~= l)
                        coefficient = coefficient  * (xi(k)-x)/(xi(k)-xi(i));
                    end  
                end
                
                derShapeFunctions(i) = derShapeFunctions(i) - 1/(xi(l)-xi(i)) * coefficient;
                
            end
        end
        
    end

end

    