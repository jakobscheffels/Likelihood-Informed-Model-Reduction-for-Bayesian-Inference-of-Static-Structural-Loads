function [posCov, G, df] = solveReducedModel(V,W)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to solve for posterior covariance using model reduction
    % IN:   V: reconstruction basis 
    %       W: projection basis
    % OUT:  posCov: posterior covariance approximation
    %       G: approximation of forward operator
    %       df: foerstner distance to true posterior covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat K B
    
    K_red = W'*K*V;
    
    G = B*V*(K_red\W');
    
    [posCov, df] = posteriorCovarianceReduced(G);

end