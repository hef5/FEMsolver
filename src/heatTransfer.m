
% dimension definition
nDimensions = size(id, 2);
nNodes = size(id, 1);
nElements = size(it, 2);
nElementNodes = size(it, 1);
nEquations  = size(id(id>0), 1);
nContraints = size(id(id<0), 1); % constraints due to dirichlet boundary condition


% global matrix initialization
boundaryLoadVector = zeros(nEquations, 1);
boundaryLoadConstrait = zeros(nContraints, 1); 
stiffnessMatrix = zeros(nEquations, nEquations);
constraitMatrix = zeros(nEquations, nContraints);

for e = 1: nElements  

    % array[nNodes][1] mesh connectivity of e-th element;
    eIt = it(:, e);
    % array[nNodes][nDimensions] equation id of e-th element;
    eId = id(eIt, :);

    for  g=1: nGauss
        for k=1: nGauss
            
            [shapeFunction, derShapeFunction]...
                = getShapeFunction(xi1(g), xi2(k));  
            jacobiMatrix... 
                = computeJacobiMatrix(derShapeFunction, nodeCoordinates);
            jacobiDeterminant... 
                = computeJacobiDeterminant(jacobiMatrix);
            der2PyhsicaCoord...
                =computeDer2PhysicalCood(derShapeFunction, jacobiMatrix);
            
            for iNode = 1: nNodes
                for jNode = 1: jNodes
                    bMatrixI = der2PyhsicaCoord(iNode, :);  % buid the B matrix
                    bMatrixJ = der2PyhsicaCoord(jNode, :);
                    matrixAtNodeIj = matrixAtNodeIj + bMatrixI' * materialMatrix * bMatrixJ *  jacobiDeterminant * B;

                    % assembly
                    [stiffnessMatrix, constraitMatrix] = assemblyMatrix(stiffnessMatrix, constraitMatrix, matrixAtNodeIj, iNode, jNode, eId, nDimensions);
                end
            end
            
        
        end
    end
    
    % stiffnessMatrix = assemblyMatrix( eStiffnessMatrix )
end
