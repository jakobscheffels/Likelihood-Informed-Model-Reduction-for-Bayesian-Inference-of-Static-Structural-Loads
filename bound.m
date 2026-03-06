clearvars

addpath src\
%% Parameters of the problem
rng(25);

beam_bool = false;
tunnel = false;

Parameters({'beam_bool','tunnel'},...
    {beam_bool,tunnel});
load Parameters.mat

%% Setup inverse problem
x_dofs=applyBoundaryCondition(x_dofs,BC_dofs,'Coordinate');

if tunnel
    LIP_Setup({'gamma_obs'},{0.05^2});
else
    LIP_Setup({'gamma_obs'},{0.001^2});
end

load('LIP_Setup.mat')
S_obs = chol(gamma_obs,'lower');
%% Generalized eigenvectors / LIS basis
[Omega,~,Nu]=svd((S_obs\G)*S_pr);
%V = S_pr*omega;
[delta,V,W,~] = calculateLISBasis();
%[~,state_samples]=gen_samples(100);

%% Calculate POD basis using 50 samples
%[~,~,Phi,~] = gen_samples(10);
%d_f_POD10=zeros(1,m);

d_f_LI=zeros(1,m);
d_f_OLR=zeros(1,m);
d_f_Sta=zeros(1,m);
p_LI=zeros(1,m);
p_OLR=zeros(1,m);
%p_POD=zeros(1,m);

R_bound=zeros(1,m);
R_bound_new=zeros(1,m);
R_boundSep=zeros(1,m);
LIS_bound = zeros(1,m);
LIS_boundSep = zeros(1,m);

% Constant
C_bound = zeros(1,m);
C_boundOLR = zeros(1,m);
HessErr = zeros(1,m);
HessErrOLR = zeros(1,m);
OLR_bound = eps.*ones(1,m);

for i=1:m

    %% LIS reduced operator
    [gamma_pos_LI,G_LI_cell{i},d_f_LI(i)]=solveReducedModel(V(:,1:i),W(:,1:i));
    
    p_LI(i)=norm(gamma_pos-gamma_pos_LI,2);

    
    
    %% Spantini Reduction
    [gamma_pos_OLR,G_OLR_cell{i},d_f_OLR(i)]=solveOLRA(V(:,1:i),W(:,1:i));
    p_OLR(i)=norm(gamma_pos-gamma_pos_OLR,2);
 

    %% LIS bound
    R = Omega(:,1:i)'*(S_obs'\C)*S_pr*Nu(:,1:i);
    R_bound(i)=norm((Omega(:,1:i)-(S_obs\C)*S_pr*(Nu(:,1:i)/R))*diag(delta(1:i)),2);
    R_boundSep(i)=norm((Omega(:,1:i)-(S_obs\C)*S_pr*(Nu(:,1:i)/R)),2);
    R_bound_new(i)=norm(Omega(:,i+1:end)'*(S_obs\C)*S_pr*Nu(:,1:i)*(R\(diag(delta(1:i)))),2);
    HessErrOLR(i)=norm((S_obs\(G-G_OLR_cell{i}))*S_pr,2);


    G_hat=C*V(:,1:i)*((W(:,1:i)'*K*V(:,1:i))\W(:,1:i)');
    HessErr(i)=norm((S_obs\(G-G_hat))*S_pr,2);
    
    preHess = (S_obs\G)*S_pr;
    preHessAp = (S_obs\G_hat)*S_pr;
    preHessOLR = (S_obs\G_hat)*S_pr;
    norm_hess = norm(preHess,2);
    norm_hess_ap = norm(preHessAp,2);
    norm_hess_OLR = norm(preHessOLR,2);

    C_bound(i) = norm(S_pr,2)^2*(norm((eye(m)+preHess*preHess')\preHess,2)+ ...
    norm_hess_ap*norm_hess*(norm_hess_ap+norm_hess) ...
    +norm((eye(m)+preHessAp*preHessAp')\preHessAp,2));
    C_boundOLR(i) = norm(S_pr,2)^2*(norm((eye(m)+preHess*preHess')\preHess,2)+ ...
    norm_hess_OLR*norm_hess*(norm_hess_OLR+norm_hess) ...
    +norm((eye(m)+preHessOLR*preHessOLR')\preHessOLR,2));

    if i< m
        LIS_bound(i)=delta(i+1)+R_bound(i);
        LIS_boundSep(i)=delta(i+1)+R_boundSep(i)*delta(1);

        OLR_bound(i)=delta(i+1);
    else
        LIS_bound(i)=eps;
        LIS_boundSep(i)=eps;
    end

end


figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
semilogy(LIS_bound,"LineWidth",2)
set(gca,"FontSize",20)
hold on
semilogy(LIS_boundSep,"--","LineWidth",2)
semilogy(HessErr,"LineWidth",2)
semilogy(OLR_bound,'x',"MarkerSize",6,"LineWidth",1.5)
semilogy(HessErrOLR,"o","MarkerSize",6,"LineWidth",1.5)
xlabel("r","Interpreter","latex")
title("$\Vert S_{obs}^{-1}(G-\widehat{G})S \Vert_p\le\delta_{r+1}+\Vert(\Omega_r-S_{obs}^{-1}CS\bar{\nu}_rR_r^{-1})\Delta_r\Vert_p$","Interpreter","latex")
axis([1 10 eps 1e2])
legend("Joint","Separate","Actual","OLR Bound","OLR error","Interpreter","latex","Location","southwest")

nexttile;
semilogy(C_bound.*LIS_bound,"LineWidth",2)
set(gca,"FontSize",20)
hold on
semilogy(C_bound.*LIS_boundSep,"--","LineWidth",2)
semilogy(p_LI,"LineWidth",2)
semilogy(C_boundOLR.*OLR_bound,"x","MarkerSize",6,"LineWidth",1.5)
semilogy(p_OLR,"o","MarkerSize",6,"LineWidth",1.5)
xlabel("r","Interpreter","latex")
title("Posterior covariance error for Bar problem","Interpreter","latex")
legend("Joint","Separate","Actual","OLR Bound","OLR error","Location","southwest","Interpreter","Latex")
axis([1 10 1e-5 1e15])


figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
semilogy(delta,"LineWidth",2)
set(gca,"FontSize",20)
title("Singular values","Interpreter","latex")
xlabel("r","Interpreter","latex")
axis([1 10 eps 1e2])

nexttile;
semilogy(R_bound,"LineWidth",2)
set(gca,"FontSize",20)
hold on
semilogy(R_boundSep.*delta(1),"--","LineWidth",2)
semilogy(R_bound_new,"x","MarkerSize",6,"LineWidth",1.5)
title("$\Vert (\Omega_r-S_{obs}^{-1}CS\bar{\nu}_rR_r^{-1})\Delta_r\Vert_p$","Interpreter","latex")
axis([1 10 eps 1e2])
xlabel("r","Interpreter","latex")
legend("Joint","Separate","New","Interpreter","latex","Location","southwest")