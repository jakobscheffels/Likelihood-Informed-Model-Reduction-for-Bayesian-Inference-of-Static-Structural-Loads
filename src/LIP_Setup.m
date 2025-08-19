function LIP_Setup(inputName, inputVal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Input: inputName is cell containing variable names to be modified
    %%%            inputVal is list containing values of variables
    %%%
    %%%     Output: Stores modified parameters in Parameter.mat to be loaded
    %%%
    %%%     Valid variable names (default value):
    %%%       nm: number of measurements (10)
    %%%       gamma_obs: variance of observation noise (1e-5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng(41)

    load('Parameters.mat')
    
    nm = 10;
    gamma_obs = 1e-5;


    if nargin > 0
        for i=1:size(inputName,2)
            input = inputVal(i);
            if strcmp(inputName(i),'nm')
                nm=input{1};
            elseif strcmp(inputName(i),'gamma_obs')
                gamma_obs=input{1};
            else
                inputStr = inputName(i);
                inputStr = inputStr{:};
                text = append('Invalid variable name: ',inputStr,'! Using default parameter instead!');
                warning(text)
            end
        end
    end
    
    nm_pos = randi([1,nnode],1,nm);
    nm_pos = unique(nm_pos);
    nm_pos = setdiff(nm_pos,BC_dofs);
    
    while ~(size(unique(nm_pos),2)==nm)
        nm_pos=[nm_pos, randi([2,nele],1,nm-size(nm_pos,2))];
        nm_pos = unique(nm_pos);
        nm_pos = setdiff(nm_pos,BC_dofs);
    end
    
    nm_pos=sort(nm_pos);
    %nm_pos(end+1)=nnode;
    %nm_pos = nm:(nnode-1)/nm:nnode;
    %nm_pos(end+1)=nnode;
    gamma_obs=gamma_obs.*diag(ones(1,nm));
    
    if beam_bool
        
        if tunnel
            EI = pi/64*zeta*E*(D^4-(D-2*t)^4);
            K = assembleBeamK(EI,1,l,nele,nnode);
            KG = assembleGroundStiffnessK(k_vector);
            K = K + KG;
        else
            K = assembleBeamK(E,I,L,nele,nnode);
        end
    else
        Ke = D/l;
        K = assembleK(Ke,nele,nele+1);
    end
    K = applyBoundaryCondition(K,BC_dofs);
    
    B = observationMapper(nm_pos);
    B = applyBoundaryCondition(B,BC_dofs,'Observation');
  
    G = (K\(B'))';

    xGL=0;
    XQ = l/2*(repmat(xGL',1,nele)+repmat(1:2:2*nele,length(xGL),1));
    xq = XQ(:);
    %phi_at_xq = cell2mat(cellfun(@(c) c(xq'),phi,'un',0));
    gamma_prior_q=zeros(nele,nele);
    for i = 1:nele
        for j=1:nele
            gamma_prior_q(i,j)=sigma_q^2*exp(-abs(xq(i)-xq(j))/theta);
        end
    end

    gamma_prior_f=L_mat*gamma_prior_q*L_mat';
    gamma_prior_f=applyBoundaryCondition(gamma_prior_f,BC_dofs);
   
    singular_prior=true;

    gamma_prior_f_inv = zeros(size(gamma_prior_f));
    S_pr = L_mat*chol(gamma_prior_q,'lower');
    S_pr = applyBoundaryCondition(S_pr,BC_dofs,'Force');

    mu_f = L_mat*(mu_q.*ones(nele,1));
    mu_f = applyBoundaryCondition(mu_f,BC_dofs,'Force');

    tmp = G*gamma_prior_f*G'+gamma_obs;
    gamma_pos = gamma_prior_f-gamma_prior_f*G'*(tmp\G)*gamma_prior_f;

    save('LIP_Setup','B','G','gamma_obs','gamma_prior_f', 'nm_pos',...
        'gamma_prior_f_inv','K','L_mat','nm','S_pr','mu_f','singular_prior','gamma_pos','gamma_prior_q')
end