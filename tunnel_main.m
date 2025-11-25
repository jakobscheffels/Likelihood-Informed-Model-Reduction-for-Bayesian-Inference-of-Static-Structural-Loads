clearvars
addpath src\
%% Parameters
rng(41)

run_POD_analysis = false; % set 'true' if want to run POD analysis

beam_bool = true;
tunnel = true;
fixed={'true'};
L = 200;
E = 35e9;
D = 6.2;
nele = 800;

BC_dofs=[];
theta = L/2;
mu_q = 300e5;
delta_q = 1;
sigma_q = mu_q*delta_q;

k_vector = zeros(nele+1,1);
k1=33000e3;
k2=5000e3;
k_vector(1:nele/2)=k1.*ones(nele/2,1);
k_vector(nele/2+1:end)=k2.*ones(nele/2+1,1);

t = 0.35;
zeta = 1/7;

Parameters({'beam_bool','BC_dofs','fixed','L','nele','E','D','theta','sigma_q','mu_q','k_vector','zeta','t','tunnel'},...
    {beam_bool,BC_dofs,fixed,L,nele,E,D,theta,sigma_q,mu_q,k_vector,zeta,t,tunnel});
load Parameters.mat

N_rep = 200;
%% Setup inverse problem

LIP_Setup();
load LIP_Setup.mat C gamma_prior_f gamma_obs gamma_pos G m K mu_f S_pr 
[~,~,omega]=svd((chol(gamma_obs,'lower')\G)*S_pr);
V = S_pr*omega;
%% POD Analysis
if run_POD_analysis
    [S_f10,S_u10,Phi10,Sigma10] = gen_samples(10);
    [S_f50,S_u50,Phi50,Sigma50] = gen_samples(50);
    [S_f100,S_u100,Phi100,Sigma100] = gen_samples(100);
    [S_f200,S_u200,Phi200,Sigma200] = gen_samples(200);
    [S_f500,S_u500,Phi500,Sigma500] = gen_samples(500);
    S_f10 = L_mat*S_f10;
    S_f50 = L_mat*S_f50;
    S_f100 = L_mat*S_f100;
    S_f200 = L_mat*S_f200;
    S_f500 = L_mat*S_f500;
    
    [~,S_uTest,~,~] = gen_samples(200);
    
    proj_errorPOD10 = zeros(1,10);
    proj_errorPOD50 = zeros(1,50);
    proj_errorPOD100 = zeros(1,100);
    proj_errorPOD200 = zeros(1,200);
    proj_errorPOD500 = zeros(1,500);
    
    for i=1:500
        if i<11
            proj_errorPOD10(i) = mean(vecnorm(S_uTest-Phi10(:,1:i)*Phi10(:,1:i)'*S_uTest,2,1));
        end
        if i<51
            proj_errorPOD50(i) = mean(vecnorm(S_uTest-Phi50(:,1:i)*Phi50(:,1:i)'*S_uTest,2,1));
        end
        if i<101
            proj_errorPOD100(i) = mean(vecnorm(S_uTest-Phi100(:,1:i)*Phi100(:,1:i)'*S_uTest,2,1));
        end
        if i<201
            proj_errorPOD200(i) = mean(vecnorm(S_uTest-Phi200(:,1:i)*Phi200(:,1:i)'*S_uTest,2,1));
        end
        proj_errorPOD500(i) = mean(vecnorm(S_uTest-Phi500(:,1:i)*Phi500(:,1:i)'*S_uTest,2,1));
    end
else
    [~,~,Phi10,~] = gen_samples(10);
end
%% Model Reduction
[delta,V,W]=calculateLISBasis();
[~,state_samples] = gen_samples(100);
save LIS_Basis_Tunnel.mat V W K mu_f state_samples
Phi = Phi10;
d_f_OLR = zeros(1,m);
d_f_LI = zeros(1,m);
d_f_POD = zeros(1,m);

for i =1:m
    [posCovLI,G_LI_cell{i},d_f_LI(i)] = solveReducedModel(V(:,1:i),W(:,1:i));

    [posCovOLR,G_OLR_cell{i},d_f_OLR(i)] = solveOLRA(V(:,1:i),W(:,1:i));

    [posCovPOD,G_POD_cell{i},d_f_POD(i)] = solveReducedModel(Phi(:,1:i),Phi(:,1:i));
end

%% Mean 
[~,state_sample] = gen_samples(N_rep);

error_LI=zeros(N_rep,m);
error_POD=zeros(N_rep,m);
error_OLR=zeros(N_rep,m);
for j=1:N_rep
    ysam = C*state_sample(:,j)+sqrt(gamma_obs)*randn(m,1);
   
    mu_full = meanCalculation(G,ysam);
    mu_full_norm = norm(mu_full);

    for i = 1:m
        mu_LI = meanCalculation(G_LI_cell{i},ysam);

        mu_POD = meanCalculation(G_POD_cell{i},ysam);

        mu_OLR = meanCalculation(G_OLR_cell{i},ysam);

        error_LI(j,i) = norm(mu_full-mu_LI)/mu_full_norm;
        error_POD(j,i) = norm(mu_full-mu_POD)/mu_full_norm;
        error_OLR(j,i) = norm(mu_full-mu_OLR)/mu_full_norm;
    end
    
end
mean_LI = mean(error_LI,1);
mean_POD = mean(error_POD,1);
mean_OLR = mean(error_OLR,1);


%% PLOT
width = 14.5;
height = 6;
if run_POD_analysis

    figure
    t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    loglog(diag(Sigma10),'--o','LineWidth',2,'MarkerSize',get(0,'DefaultLineMarkerSize'))
    box off
    set(gca,'FontSize',20)
    hold on
    loglog(diag(Sigma50),'LineWidth',2,'LineStyle','--')
    loglog(diag(Sigma100),'LineWidth',2,'LineStyle',':')
    loglog(diag(Sigma200),'LineWidth',2,'LineStyle','-.')
    loglog(diag(Sigma500),'LineWidth',2)
    xlabel('Approximation rank r', 'Interpreter','latex')
    title('Singular value decay','Interpreter','latex','FontSize',28)
    legend('10','50','100','200','500','Location','southwest')
    legend boxoff
    axis([1 10^3 10^-10 10^4])
 
    nexttile;
    
    loglog(proj_errorPOD10,'--o','LineWidth',2,'MarkerSize',get(0,'DefaultLineMarkerSize'))
    box off
    set(gca,'FontSize',20)
    hold on
    loglog(proj_errorPOD50,'LineWidth',2,'LineStyle','--')
    loglog(proj_errorPOD100,'LineWidth',2,'LineStyle',':')
    loglog(proj_errorPOD200,'LineWidth',2,'LineStyle','-.')
    loglog(proj_errorPOD500,'LineWidth',2)
    legend('10','50','100','200','500','Location','southwest')
    legend boxoff
    title('Projection error of POD','Interpreter','latex','FontSize',28)
    xlabel('Approximation rank r', 'Interpreter','latex')
    axis([1 10^3 10^-10 10^4])
    
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [.5 .5 width height]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [width height]);
    set(gcf, 'PaperPosition', [0 0 width height]);

    exportgraphics(gcf, 'PODTunnel.pdf', 'ContentType', 'vector');
    
end

height= 8;
alpha = 0.5;
LI_color = (1-alpha)*[0.4660 0.6740 0.1880]+alpha*[1 1 1];
OLR_color = (1-alpha)*[0.8500 0.3250 0.0980]+alpha*[1 1 1];
alpha=0.25;
POD_color = (1-alpha)*[0.3010 0.7450 0.9330]+alpha*[1 1 1];

figure
t = tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

% First plot
ax1 = nexttile;
imagesc(x_dofs,x_dofs,gamma_prior_f(1:2:end,1:2:end))
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ylabel('z','Interpreter','latex')
title("Prior covariance",'Interpreter','latex','FontSize',28)

ax2 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos(1:2:end,1:2:end))
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f(1:2:end,1:2:end)))];
yticks([])
%title('Analytical posterior covariance $\mathbf{\Gamma}_{\mathrm{pos}}$','Interpreter','latex','FontSize',18)
title('Posterior covariance','Interpreter','latex','FontSize',28)

ax3 = nexttile;
imagesc(x_dofs,x_dofs,posCovLI(1:2:end,1:2:end))
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f(1:2:end,1:2:end)))];
yticks([])
colorbar
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\mathrm{LIS}}$','Interpreter','latex','FontSize',28)

% Add colorbar below both plots using the first axes
cb = colorbar(ax3, 'Location', 'eastoutside');

nexttile;
axis off;

ax5 = nexttile;
imagesc(x_dofs,x_dofs,posCovOLR(1:2:end,1:2:end))
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ylabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
%yticks([])
%colorbar
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\mathrm{OLR}}$','Interpreter','latex','FontSize',28)

ax6 = nexttile;
imagesc(x_dofs,x_dofs,posCovPOD(1:2:end,1:2:end))
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
colorbar
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\mathrm{POD}}$','Interpreter','latex','FontSize',28)


set(gcf, 'Units', 'inches');
set(gcf, 'Position', [.5 .5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

exportgraphics(gcf, 'priorTunnel.pdf', 'ContentType', 'vector');

height = 6;

figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
semilogy(mean_LI,'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',20)
box off
hold on
semilogy(mean_POD,'Color',POD_color,'LineWidth',2)
semilogy(mean_OLR,'o','Color',OLR_color,'LineWidth',2)

legend('LIS','POD','OLR','Location','southwest')
legend boxoff
title('Relative posterior mean error','Interpreter','latex','FontSize',28)
axis([1 10 1e-15 100])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])
xlabel('Approximation rank $r$','Interpreter','latex')

% Second plot
ax = nexttile;
semilogy(sqrt(d_f_LI),'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',20)
box off
hold on
semilogy(sqrt(d_f_POD),'Color',POD_color,'LineWidth',2)
semilogy(sqrt(d_f_OLR),'o','Color',OLR_color,'LineWidth',2)
title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',28)
%ylabel('F$\ddot{o}$rstner distance','Interpreter','latex')
xlabel('Approximation rank $r$','Interpreter','latex')
legend('LIS','POD','OLR','Location','southwest')
legend boxoff
axis([1 10 1e-15 100])
yticks([10^(-15) 10^(-10) 10^(-5) 1])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [.5 .5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

exportgraphics(gcf, 'posTunnel.pdf', 'ContentType', 'vector');

theta=zeros(1,10);
for i=1:10
    theta(i)=1/(norm(V(:,i))*norm(W(:,i)));
end