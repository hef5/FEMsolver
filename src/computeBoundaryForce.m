

function Re = computeBoundaryForce(p, elementOfSideSet, sideOfSideSet, x_e, EQ_e, nElements, t, it)

nElementNodes = size(it, 1);
% elemOfSideSet: list of the element on neumann boundary.
% sideOfSideSet:
% nodeOfSideSet:

Re = zeros(nElementNodes, 2, nElements);

xi = gaussPoint(p);
alpha = gaussWeight(p);

% y : NG x NN (number of gauss points x number of nodes)
% dy: NG x NN

% t:  DF_e x ND x NE_b
[t_e] = transf_e (t, it, nElementNodes, nElements);
% t_e : ND x NN x NE
NG = length(xi);

% numLocalEquation = reshape(1:EQ_e, 2, []);

for e = 1 : length(elementOfSideSet)
    
    nodeOfSideSet = getNodeSideSet (sideOfSideSet(e), p);
    elCoordinates = it(nodeOfSideSet, e);
    %numNodeGloble = it(nodeOfSideSet, elemOfSideSet(e));

    for g = 1: NG
        [shapeFunctions, derShapeFunctions] = evaluateShapeFunctions1D(p, xi(g));
        
        % 2d problem, we need line integration on the Neumman boundary,
        % therefore, 1D shape functions are used and there are only two
        % component in Jacobi matrix
        jacobiMatrix = derShapeFunctions(g, :) * elCoordinates;
        jacobiDeterminant = sqrt(jacobiMatrix(1).^2+jacobiMatrix(2).^2);
        weight = alpha(g) * jacobiDeterminant;
        
        % computing the load on the boundary through load functions
        t    = shapeFunctions(g, :) * t_e( :, nodeOfSideSet, elementOfSideSet(e))';

        Re(nodeOfSideSet, :, elementOfSideSet(e)) = Re(nodeOfSideSet, :, elementOfSideSet(e)) + shapeFunctions(g, :)' * t * weight;
    end
    
    
    
end