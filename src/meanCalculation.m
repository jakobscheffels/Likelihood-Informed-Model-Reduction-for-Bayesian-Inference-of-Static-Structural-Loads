function mu = meanCalculation(G,ysam)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate posterior mean
    % IN:   G: forward operator of inverse problem
    %       ysam: measurement sample
    % OUT:  mu: posterior mean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat 'mu_f' 'gamma_prior_f' 'gamma_obs'

    mu = mu_f+gamma_prior_f*G'*((G*gamma_prior_f*G'+gamma_obs)\(ysam-G*mu_f));
end