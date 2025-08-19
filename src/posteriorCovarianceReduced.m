function [posCov, df] = posteriorCovarianceReduced(G)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to evaluate posterior covariance and foersnter distance
    % IN:   G: forward operator of inverse problem
    %       
    % OUT:  posCov: posterior covariance for given inverse problem
    %       df: foerstner distance to true posterior covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat gamma_obs gamma_prior_f_inv singular_prior gamma_prior_f
    load LIS_Basis.mat V W

    prior_hat = gamma_prior_f;
    posCov = prior_hat-prior_hat*G'*((G*prior_hat*G'+gamma_obs)\G)*prior_hat;

    df = foerstnerDistance(posCov);
end

function df = foerstnerDistance(gamma1)
    load LIP_Setup.mat gamma_pos
    load LIS_Basis.mat W
  
  
    L = chol(W'*gamma1*W,'lower');
    R = chol(W'*gamma_pos*W,'lower');
    
    [df,~]=gen_eigenvalue(L,R);
    df(abs(df)<eps)=[];
    df = dot(log(df),log(df));
end

function [delta, w] = gen_eigenvalue(L,R)

    tmp = L'/R';
    [U,S,~]=svd(tmp);
    delta = diag(S).^2;
    w = R*U;

end