load LIS_Basis_Tunnel.mat 
V_tun = V;
W_tun = W;
K_tun = K;
mu_f_tun = mu_f;
N_tun = state_samples;
force_tun = force_samples;
C_tun = C;
G_tun = G;
d_p_tun = d_p;
gamma_prior_f_tun = gamma_prior_f;
S_pr_tun = S_pr;
Phi_tun = Phi;

load LIS_Basis_Bar.mat
V_bar = V;
W_bar = W;
K_bar = K;
mu_f_bar = mu_f;
N_bar = state_samples;
C_bar = C;
G_bar = G;
S_pr_bar = S_pr;
gamma_prior_f_bar = gamma_prior_f;
d_p_bar = d_p;
force_bar = force_samples;
Phi_bar = Phi;

load LIS_Basis_Beam.mat
V_beam = V;
W_beam=W;
K_beam=K;
mu_f_beam=mu_f;
N_beam = state_samples;
C_beam = C;

S_obs=sqrt(gamma_obs);

sigma_max_bar = zeros(1,10);
sigma_max_beam = zeros(1,10);
sigma_max_tun = zeros(1,10);
sigma_max_hess = zeros(1,10);
for i= 1:10
    [~,S,~]=svd(V_bar(:,1:i)*W_bar(:,1:i)');
    sigma_max_bar(i)=S(1,1);

    [~,S,~]=svd(V_beam(:,1:i)*W_beam(:,1:i)');
    sigma_max_beam(i)=S(1,1);

    [~,S,~]=svd(V_tun(:,1:i)*W_tun(:,1:i)');
    sigma_max_tun(i)=S(1,1);
end

P_v_bar = zeros(size(V_bar,1),size(V_bar,1));
P_v_beam = zeros(size(V_beam,1),size(V_beam,1));
P_v_tun = zeros(size(V_tun,1),size(V_tun,1));
si_bar = zeros(size(N_bar,2),10);
si_beam = zeros(size(N_beam,2),10);
si_tun = zeros(size(N_tun,2),10);
u_bar = N_bar;
u_beam = N_beam;
u_tun = N_tun;

si_bar_pod = zeros(size(N_bar,2),10);
si_tun_pod = zeros(size(N_tun,2),10);
%u_bar = (K_bar\mu_f_bar);
%u_beam = (K_beam\mu_f_beam);
%u_tun = K_tun\mu_f_tun;
c = zeros(200,10);

for i=1:10
    P_v_bar=V_bar(:,1:i)*((V_bar(:,1:i)'*V_bar(:,1:i))\V_bar(:,1:i)');
    P_v_beam=V_beam(:,1:i)*((V_beam(:,1:i)'*V_beam(:,1:i))\V_beam(:,1:i)');
    P_v_tun=V_tun(:,1:i)*((V_tun(:,1:i)'*V_tun(:,1:i))\V_tun(:,1:i)');

    P_v_bar_pod = Phi_bar(:,1:i)*Phi_bar(:,1:i)';
    P_v_tun_pod = Phi_tun(:,1:i)*Phi_tun(:,1:i)';
    c(:,i) = P_v_beam*N_beam(:,1);
    for j=1:size(N_bar,2)

        co = (N_bar(:,j)'*P_v_bar*N_bar(:,j))/(norm(N_bar(:,j))*norm(P_v_bar*N_bar(:,j)));
        si_bar(j,i) = sqrt(1-co^2);
    
        
        co = (N_beam(:,j)'*P_v_beam*N_beam(:,j))/(norm(N_beam(:,j))*norm(P_v_beam*N_beam(:,j)));
        si_beam(j,i) = sqrt(1-co^2);
    
        
        co = (N_tun(:,j)'*P_v_tun*N_tun(:,j))/(norm(N_tun(:,j))*norm(P_v_tun*N_tun(:,j)));
        si_tun(j,i) = sqrt(1-co^2);

        %POD 
        co = (N_bar(:,j)'*P_v_bar_pod*N_bar(:,j))/(norm(N_bar(:,j))*norm(P_v_bar_pod*N_bar(:,j)));
        si_bar_pod(j,i) = sqrt(1-co^2);
        co = (N_tun(:,j)'*P_v_tun_pod*N_tun(:,j))/(norm(N_tun(:,j))*norm(P_v_tun_pod*N_tun(:,j)));
        si_tun_pod(j,i) = sqrt(1-co^2);
    end
    %for j=101:400
    %    co = (N_tun(:,j)'*P_v_tun*N_tun(:,j))/(norm(N_tun(:,j))*norm(P_v_tun*N_tun(:,j)));
    %    si_tun(j,i) = sqrt(1-co^2);
    %end
end

%% Calculate E[|u-u_r|^2] and E[|C(u-u_r)|^2]
u_diff = zeros(1,10);
y_diff = zeros(1,10);
y_diff_tun = zeros(1,10);
u_diff_pod = zeros(1,10);
u_diff_tun = zeros(1,10);
y_diff_pod = zeros(1,10);
y_diff_tun_pod = zeros(1,10);
u_diff_tun_pod = zeros(1,10);
lower_bound_bar = zeros(1,10);
upper_bound_bar = zeros(1,10);
y_lower_bound_bar = zeros(1,10);
y_upper_bound_bar = zeros(1,10);
lower_bound_bar_pod = zeros(1,10);
upper_bound_bar_pod = zeros(1,10);
y_lower_bound_bar_pod = zeros(1,10);
y_upper_bound_bar_pod = zeros(1,10);

lower_bound_tun = zeros(1,10);
upper_bound_tun = zeros(1,10);
y_lower_bound_tun = zeros(1,10);
y_upper_bound_tun = zeros(1,10);
lower_bound_tun_pod = zeros(1,10);
upper_bound_tun_pod = zeros(1,10);
y_lower_bound_tun_pod = zeros(1,10);
y_upper_bound_tun_pod = zeros(1,10);

alpha=0.5;

for j=1:10
    v = V_bar(:,1:j);
    w = W_bar(:,1:j);
    K_red = w'*K_bar*v;
    u_r = v*(K_red\(w'*force_bar));
    e = vecnorm(N_bar-u_r,2,1).^2;
    y = vecnorm(C_bar*(N_bar-u_r),2,1).^2;
    u_diff(j)=mean(e);
    y_diff(j)=mean(y);

    s_e = std(e)/sqrt(size(N_bar,2));

    tval = tinv(1-alpha/2,size(N_bar,2)-1);

    lower_bound_bar(j)=u_diff(j)-tval*s_e;
    upper_bound_bar(j)=u_diff(j)+tval*s_e;

    s_y = std(y)/sqrt(size(N_bar,2));

    %tval = tinv(1-alpha/2,size(N_bar,2)-1);

    y_lower_bound_bar(j)=y_diff(j)-tval*s_y;
    y_upper_bound_bar(j)=y_diff(j)+tval*s_y;
    
    v = V_tun(:,1:j);
    w = W_tun(:,1:j);
    K_red = w'*K_tun*v;
    u_r = v*(K_red\(w'*force_tun));
    e = vecnorm(N_tun-u_r,2,1).^2;
    y = vecnorm(C_tun*(N_tun-u_r),2,1).^2;
    u_diff_tun(j)=mean(e);
    y_diff_tun(j) = mean(y);
    s_e = std(e)/sqrt(size(N_tun,2));

    tval = tinv(1-alpha/2,size(N_tun,2)-1);

    lower_bound_tun(j)=u_diff_tun(j)-tval*s_e;
    upper_bound_tun(j)=u_diff_tun(j)+tval*s_e;

    s_y = std(y)/sqrt(size(N_tun,2));

    %tval = tinv(1-alpha/2,size(N_tun,2)-1);

    y_lower_bound_tun(j)=y_diff_tun(j)-tval*s_y;
    y_upper_bound_tun(j)=y_diff_tun(j)+tval*s_y;

    %POD
    v = Phi_bar(:,1:j);
    w = Phi_bar(:,1:j);
    K_red = w'*K_bar*v;
    u_r = v*(K_red\(w'*force_bar));
    e = vecnorm(N_bar-u_r,2,1).^2;
    y = vecnorm(C_bar*(N_bar-u_r),2,1).^2;
    u_diff_pod(j)=mean(e);
    y_diff_pod(j)=mean(y);

    s_e = std(e)/sqrt(size(N_bar,2));

    tval = tinv(1-alpha/2,size(N_bar,2)-1);

    lower_bound_bar_pod(j)=u_diff_pod(j)-tval*s_e;
    upper_bound_bar_pod(j)=u_diff_pod(j)+tval*s_e;

    s_y = std(y)/sqrt(size(N_bar,2));

    %tval = tinv(1-alpha/2,size(N_bar,2)-1);

    y_lower_bound_bar_pod(j)=y_diff_pod(j)-tval*s_y;
    y_upper_bound_bar_pod(j)=y_diff_pod(j)+tval*s_y;
    
    v = Phi_tun(:,1:j);
    w = Phi_tun(:,1:j);
    K_red = w'*K_tun*v;
    u_r = v*(K_red\(w'*force_tun));
    e = vecnorm(N_tun-u_r,2,1).^2;
    y = vecnorm(C_tun*(N_tun-u_r),2,1).^2;
    u_diff_tun_pod(j)=mean(e);
    y_diff_tun_pod(j)=mean(y);

    s_e = std(e)/sqrt(size(N_tun,2));

    tval = tinv(1-alpha/2,size(N_tun,2)-1);

    lower_bound_tun_pod(j)=u_diff_tun_pod(j)-tval*s_e;
    upper_bound_tun_pod(j)=u_diff_tun_pod(j)+tval*s_e;

    s_y = std(y)/sqrt(size(N_tun,2));

    %tval = tinv(1-alpha/2,size(N_tun,2)-1);

    y_lower_bound_tun_pod(j)=y_diff_tun_pod(j)-tval*s_y;
    y_upper_bound_tun_pod(j)=y_diff_tun_pod(j)+tval*s_y;

end
%{
figure
subplot(1,2,1)
plot(N_beam(1:2:end,1))
hold on
plot(c(1:2:end,1))
plot(c(1:2:end,2))
plot(c(1:2:end,3))
plot(c(1:2:end,4))
plot(c(1:2:end,5))
title('Transversal DoFs')
legend('Full state','1','2','3','4','5','Location','northwest')

subplot(1,2,2)
plot(N_beam(2:2:end,1))
hold on
plot(c(2:2:end,1))
plot(c(2:2:end,2))
plot(c(2:2:end,3))
plot(c(2:2:end,4))
plot(c(2:2:end,5))
title('Rotational DoFs')
%}
[~,S,~]=svd((K_bar\eye(100)));
S_bar = S(1,1);
[~,S,~]=svd((K_beam\eye(200)));
S_beam = S(1,1);
[~,S,~]=svd((K_tun\eye(1602)));
S_tun = S(1,1);

%% Sprungk
y_bar = C_bar*u_bar(:,1);
factor = zeros(1,10);
Z = 1/sqrt((2*pi)^10*det(G*gamma_prior_f*G'+gamma_obs))*exp(-0.5*(y_bar-G_bar*mu_f_bar)'*((G*gamma_prior_f*G'+gamma_obs)\(y_bar-G_bar*mu_f_bar)));
for i=1:10
    G_hat = C_bar*V_bar(:,1:i)*((W_bar(:,1:i)'*K_bar*V_bar(:,1:i))\W_bar(:,1:i)');
    factor(i) = 2*sqrt((y_bar'*(gamma_obs\y_bar))+sqrt(trace(G_bar'*(gamma_obs\G_bar)))+sqrt(trace(G_hat'*(gamma_obs\G_hat))))*...
                sigma_max_bar(i)*si_bar(1,i)*S_bar;
end
% E[sin^2 *norm(u)^2]
bound_bar = mean((si_bar.^2).*(vecnorm(N_bar,2,1).^2'),1);
bound_bar_pod = mean((si_bar_pod.^2).*(vecnorm(N_bar,2,1).^2'),1);
bound_beam = mean((si_beam.^2).*(vecnorm(N_beam,2,1).^2'),1);
bound_tun = mean((si_tun.^2).*(vecnorm(N_tun,2,1).^2'),1);
bound_tun_pod = mean((si_tun_pod.^2).*(vecnorm(N_tun,2,1).^2'),1);
%% Han 
S_obs = sqrt(gamma_obs);
[~,S,~]=svd((eye(10)+(S_obs\G_bar)*gamma_prior_f_bar*(G_bar'/S_obs))\(S_obs\G_bar)*S_pr_bar);
S_1_bar = S(1,1);
[~,S,~]=svd((S_obs\G_bar)*S_pr_bar);
S_2_1_bar = S(1,1);
[~,S,~]=svd(S_pr_bar);
S_prior_bar = S(1,1);
[~,S,~]=svd(S_obs\eye(10));
S_obs_bar = S(1,1);
[~,S,~]=svd(K_bar\eye(100));
S_k_bar = S(1,1);

[~,S,~]=svd((eye(10)+(S_obs\G_tun)*gamma_prior_f_tun*(G_tun'/S_obs))\(S_obs\G_tun)*S_pr_tun);
S_1_tun = S(1,1);
[~,S,~]=svd((S_obs\G_tun)*S_pr_tun);
S_2_1_tun = S(1,1);
[~,S,~]=svd(S_pr_tun);
S_prior_tun = S(1,1);
[~,S,~]=svd(S_obs\eye(10));
S_obs_tun = S(1,1);
[~,S,~]=svd(K_tun\eye(1602));
S_k_tun = S(1,1);

Cov_bound_bar=zeros(1,10);
Cov_ms_bar = zeros(1,10);
Cov_bound_tun=zeros(1,10);
Cov_ms_tun = zeros(1,10);
for i = 1:10
    G_hat = C_bar*V_bar(:,1:i)*((W_bar(:,1:i)'*K_bar*V_bar(:,1:i))\W_bar(:,1:i)');
    
    [~,S,~]=svd((S_obs\G_hat)*S_pr_bar);
    S_2_2 = S(1,1);
    S_2 = S_2_1_bar*S_2_2*(S_2_2+S_2_1_bar);
    [~,S,~]=svd((eye(10)+(S_obs\G_hat)*gamma_prior_f_bar*(G_hat'/S_obs))\(S_obs\G_hat)*S_pr_bar);
    S_3 = S(1,1);
    Con_Bar = S_1_bar+S_2+S_3;
    
    Cov_bound_bar(i) = Con_Bar*S_prior_bar^2*S_obs_bar*sigma_max_bar(i)*si_bar(1,i)*S_k_bar*S_prior_bar;
    Cov_ms_bar(i) = Con_Bar^2*S_prior_bar^4*S_obs_bar^2*sigma_max_bar(i).^2*mean(si_bar(:,i).^2)*S_prior_bar^2;

    G_hat = C_tun*V_tun(:,1:i)*((W_tun(:,1:i)'*K_tun*V_tun(:,1:i))\W_tun(:,1:i)');
    
    [~,S,~]=svd((S_obs\G_hat)*S_pr_tun);
    S_2_2 = S(1,1);
    S_2 = S_2_1_tun*S_2_2*(S_2_2+S_2_1_tun);
    [~,S,~]=svd((eye(10)+(S_obs\G_hat)*gamma_prior_f_tun*(G_hat'/S_obs))\(S_obs\G_hat)*S_pr_tun);
    S_3 = S(1,1);
    Con_tun = S_1_tun+S_2+S_3;
    
    Cov_bound_tun(i) = Con_tun*S_prior_tun^2*S_obs_tun*sigma_max_tun(i)*si_tun(1,i)*S_k_tun*S_prior_tun;
    Cov_ms_tun(i) = Con_tun^2*S_prior_tun^4*S_obs_tun^2*sigma_max_tun(i).^2*mean(si_tun(:,i).^2)*S_prior_tun^2;

end
%% Plot
figure
semilogy(Cov_bound_bar)
hold on
semilogy(d_p_bar)
legend('Bound','Actual','Interpreter','latex','Location','southwest')
title('Error bound bar for one realization','Interpreter','latex')
ylabel('$E[\Vert \Gamma_{pos}-\hat{\Gamma}_{pos}\Vert_\infty]$','Interpreter','latex')

figure
semilogy(Cov_ms_bar)
hold on
semilogy(d_p_bar.^2)
legend('Bound','Actual','Interpreter','latex','Location','southwest')
title('Mean Square error bar','Interpreter','latex')
ylabel('$E[\Vert \Gamma_{pos}-\hat{\Gamma}_{pos}\Vert_\infty^2]$','Interpreter','latex')

% Tunnel
figure
semilogy(Cov_bound_tun)
hold on
semilogy(d_p_tun)
legend('Bound','Actual','Interpreter','latex','Location','southwest')
title('Error bound tunnel for one realization','Interpreter','latex')
ylabel('$E[\Vert \Gamma_{pos}-\hat{\Gamma}_{pos}\Vert_\infty]$','Interpreter','latex')

figure
semilogy(Cov_ms_tun)
hold on
semilogy(d_p_tun.^2)
legend('Bound','Actual','Interpreter','latex','Location','southwest')
title('Mean Square error tunnel','Interpreter','latex')
ylabel('$E[\Vert \Gamma_{pos}-\hat{\Gamma}_{pos}\Vert_\infty^2]$','Interpreter','latex')



figure
plot(sigma_max_bar.^2.*bound_bar)
hold on
plot(u_diff)
plot(lower_bound_bar)
plot(upper_bound_bar)
legend('bound','actual','','','Interpreter','Latex')
title('$E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')

figure
plot(bound_bar_pod)
hold on
plot(u_diff_pod)
plot(lower_bound_bar_pod)
plot(upper_bound_bar_pod)
legend('bound','actual','','','Interpreter','Latex')
title('POD: $E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')

figure
semilogy(bound_bar_pod)
hold on
semilogy(u_diff_pod)
semilogy(lower_bound_bar_pod)
semilogy(upper_bound_bar_pod)
legend('bound','actual','','','Interpreter','Latex')
title('POD: $E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')

figure
subplot(1,2,1)
semilogy(sigma_max_bar.^2.*bound_bar)
hold on
semilogy(u_diff)
semilogy(lower_bound_bar)
semilogy(upper_bound_bar)
legend('bound','actual','','','Interpreter','Latex')
title('$E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')
axis([1 10 1e-10 1e0])
subplot(1,2,2)
semilogy(bound_bar_pod)
hold on
semilogy(u_diff_pod)
semilogy(lower_bound_bar_pod)
semilogy(upper_bound_bar_pod)
legend('bound','actual','','','Interpreter','Latex')
title('POD: $E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')
axis([1 10 1e-10 1e0])

figure
subplot(1,2,1)
semilogy(sigma_max_tun.^2.*bound_tun)
hold on
semilogy(u_diff_tun)
semilogy(lower_bound_tun)
semilogy(upper_bound_tun)
legend('bound','actual','','','Interpreter','Latex')
title('$E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')
axis([1 10 1e0 1e3])
subplot(1,2,2)
semilogy(bound_tun_pod)
hold on
semilogy(u_diff_tun_pod)
semilogy(lower_bound_tun_pod)
semilogy(upper_bound_tun_pod)
legend('bound','actual','','','Interpreter','Latex')
title('POD: $E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')
axis([1 10 1e0 1e3])

%% Plot of E[|y-y_r|^2]
figure
subplot(1,2,1)
semilogy(y_diff)
hold on
semilogy(y_diff_pod)
title('Bar: $E[\Vert y-y_r\Vert^2]$','Interpreter','latex')
axis([1 10 1e-30 1e1])
subplot(1,2,2)
semilogy(y_diff_tun)
hold on
semilogy(y_diff_tun_pod)
title('Tunnel: $E[\Vert y-y_r\Vert^2]$','Interpreter','latex')
axis([1 10 1e-30 1e1])
legend('LIS','POD','location','southwest')

figure
plot(sigma_max_tun.^2.*bound_tun)
hold on
plot(u_diff_tun)
plot(lower_bound_tun)
plot(upper_bound_tun)
legend('bound','actual','','','Interpreter','Latex')
title('$E[\Vert u-u_r \Vert^2]\le\sigma_{max}^2(VW^\top)E[sin(u,\mathcal{P}_\mathcal{V}u)^2\Vert u\Vert ^2]$','Interpreter','latex')
xlabel('r','Interpreter','latex')

figure
plot(mean(si_bar,1))
hold on
plot(mean(si_beam,1))
plot(mean(si_tun,1))
legend('Bar','Beam','Tunnel','Interpreter','Latex')
title('Sin$(u,\mathcal{P}_\mathcal{V}u)$ for $\mu_f$','Interpreter','latex')
xlabel('r','Interpreter','latex')


figure
plot(sigma_max_bar)
hold on
plot(sigma_max_beam)
plot(sigma_max_tun)
xlabel('r','Interpreter','latex')
ylabel('$\sigma_{max}(V_rW_r^\top)$','Interpreter','latex')
title('Value of largest singular value of oblique projector','Interpreter','latex')
legend('Bar','Beam','Tunnel','Interpreter','Latex')

figure
plot(sigma_max_bar.*mean(si_bar,1))
hold on
%plot(sigma_max_beam.*si_beam)
plot(sigma_max_tun.*mean(si_tun,1))
xlabel('r','Interpreter','latex')
ylabel('$\sigma_{max}(V_rW_r^\top)E[\Vert sin(u,\mathcal{P}_\mathcal{V}u)\Vert]$','Interpreter','latex')
title('Mean of constant of error bound','Interpreter','latex')
legend('Bar','Tunnel','Interpreter','Latex')


figure
plot(S_bar.*sigma_max_bar.*si_bar)
hold on
%plot(S_beam.*sigma_max_beam.*si_beam)
plot(S_tun.*sigma_max_tun.*si_tun)
xlabel('r','Interpreter','latex')
ylabel('$\sigma_{max}(V_rW_r^\top)sin(u,\mathcal{P}_\mathcal{V}u)\Vert u\Vert_{op}$','Interpreter','latex')
title('Error bound at $\mu_f$','Interpreter','latex')
legend('Bar','Tunnel','Interpreter','Latex')

