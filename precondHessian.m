clear all

load LIS_Basis_Bar.mat

r=5;
%% State space reduction 
R = G'.*(1./sqrt(diag(gamma_obs)))';
[U,delta,Z]=svd(R'*S_pr);
delta = diag(delta);
tol = max(size(R'*S_pr))*eps(max(delta));
%r = sum(delta>tol);
U = U(:,1:r);
delta=delta(1:r);
Z = Z(:,1:r);
C_inv = C'/(C*C');

V_State = C_inv*sqrt(gamma_obs)*U;
K_hat_sta = W(:,1:r)'*K*V_State;
G_state = C*V_State*(K_hat_sta\W(:,1:r)');
Hess_Sta = S_pr'*G_state'*(gamma_obs\G_state)*S_pr;
[U_Sta,S_Sta,v_Sta]=svd(Hess_Sta);
S_Sta = diag(S_Sta);

%% Prior preconditioned Hessian
Hess = S_pr'*G'*(gamma_obs\G)*S_pr;
[U,S,v]=svd(Hess);
S = diag(S);

K_hat = W(:,1:r)'*K*V(:,1:r);
G_hat = C*V(:,1:r)*(K_hat\W(:,1:r)');
Hess_LIS = S_pr'*G_hat'*(gamma_obs\G_hat)*S_pr;
[U_LIS,S_LIS,v_LIS]=svd(Hess_LIS);
S_LIS = diag(S_LIS);

G_OLR = C*(K\V(:,1:r))*W(:,1:r)';
Hess_OLR = S_pr'*G_OLR'*(gamma_obs\G_OLR)*S_pr;
[U_OLR,S_OLR,v_OLR]=svd(Hess_OLR);
S_OLR = diag(S_OLR);

err_LIS = abs(S(1:r)-S_LIS(1:r));
err_OLR = abs(S(1:r)-S_OLR(1:r));
err_Sta = abs(S(1:r)-S_Sta(1:r));

figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;
semilogy(S(1:r),'LineWidth',2);
hold on
semilogy(S_LIS(1:r),'LineWidth',2);
semilogy(S_OLR(1:r),'LineWidth',2);
semilogy(S_Sta(1:r),'LineWidth',2);
legend('Full','LIS','OLR','State','Interpreter','latex','FontSize',20)
xlabel('r','Interpreter','latex','FontSize',20)
title('Singular values of $S^\top G^\top \Gamma_{obs}^{-1} GS$ for $r=5$','Interpreter','latex','FontSize',24)
axis([1 10 1e-18 1e2])

nexttile;
semilogy(err_LIS,'LineWidth',2)
hold on
semilogy(err_OLR,'LineWidth',2)
semilogy(err_Sta,'LineWidth',2)
legend('LIS','OLR','State','Interpreter','latex','FontSize',20)
title('Error in singular values $\vert \delta_i-\hat{\delta}_i\vert$ for $r=5$','Interpreter','latex','FontSize',24)
xlabel('r','Interpreter','latex','FontSize',20)
axis([1 10 1e-18 1e2])

figure
tiledlayout(1,2,"TileSpacing","compact",Padding="compact");
nexttile;
imagesc(K_hat)
set(gca,'FontSize',20)
title('Structure of $\hat{K}_{LIS}$','Interpreter','latex','FontSize',24)
colorbar

nexttile;
imagesc(K_hat_sta)
set(gca,'FontSize',20)
colorbar
title('Structure of $\hat{K}_{State}$','Interpreter','latex','FontSize',24)