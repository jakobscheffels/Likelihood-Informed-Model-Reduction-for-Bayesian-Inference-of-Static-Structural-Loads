clearvars

addpath src\
%% Parameters of the problem
rng(25);
run_POD_analysis = false; % set 'true' if want to run POD analysis

%% Setup of the problem
Parameters({'beam_bool'},{false});
load('Parameters.mat')

x_dofs=applyBoundaryCondition(x_dofs,BC_dofs,'Coordinate');

LIP_Setup({'gamma_obs'},{0.001^2});
load('LIP_Setup.mat')

%% Generalized eigenvectors / LIS basis
%[~,~,omega]=svd((chol(gamma_obs,'lower')\G)*S_pr);
%V = S_pr*omega;
[~,V,W,V_State] = calculateLISBasis();
%[~,state_samples]=gen_samples(100);
%save LIS_Basis_Beam.mat V W K mu_f state_samples

%% Mean samples
N_rep = 200;

% Generate sample
[~,state_sample]=gen_samples(N_rep);

if run_POD_analysis
    d_f_POD10 = zeros(1,m);
    d_f_POD20 = zeros(1,m);
    d_f_POD50 = zeros(1,m);
    d_f_POD1000 = zeros(1,m);

    %% Calculate POD basis
    
    [~,~,Phi10,~] = gen_samples(10);
    [~,~,Phi20,~] = gen_samples(20);
    [~,~,Phi50,~] = gen_samples(50);
    [~,~,Phi1000,~] = gen_samples(1000);

    for i=1:m
        [~,G_POD_cell50{i},d_f_POD50(i)]=solveReducedModel(Phi50(:,1:i),Phi50(:,1:i));

        [gamma_pos_POD,G_POD_cell10{i},d_f_POD10(i)]=solveReducedModel(Phi10(:,1:i),Phi10(:,1:i));
        [~,G_POD_cell20{i},d_f_POD20(i)]=solveReducedModel(Phi20(:,1:i),Phi20(:,1:i));
        [~,G_POD_cell1000{i},d_f_POD1000(i)]=solveReducedModel(Phi1000(:,1:i),Phi1000(:,1:i));
  
    end

    error_POD10 = zeros(N_rep,m);
    error_POD20 = zeros(N_rep,m);
    error_POD50 = zeros(N_rep,m);
    error_POD1000 = zeros(N_rep,m);

    for j=1:N_rep
        ysam = C*state_sample(:,j)+sqrt(gamma_obs)*randn(m,1);
        % Full model mean
        mu_full = meanCalculation(G,ysam);
        mu_full_norm = norm(mu_full);
    
        for i = 1:m
    
            % POD approximation
            mu_POD50 = meanCalculationPOD(G_POD_cell50{i},ysam,Phi50(:,1:i));
            mu_POD10 = meanCalculationPOD(G_POD_cell10{i},ysam,Phi10(:,1:i));
            mu_POD20 = meanCalculationPOD(G_POD_cell20{i},ysam,Phi20(:,1:i));
            mu_POD1000 = meanCalculationPOD(G_POD_cell1000{i},ysam,Phi1000(:,1:i));
  
  
            error_POD50(j,i) = norm(mu_full-mu_POD50)/mu_full_norm;
    
            error_POD10(j,i) = norm(mu_full-mu_POD10)/mu_full_norm;
            error_POD20(j,i) = norm(mu_full-mu_POD20)/mu_full_norm;
            error_POD1000(j,i) = norm(mu_full-mu_POD1000)/mu_full_norm;
        end
        
    end    
    mean_POD10 = mean(error_POD10,1);
    mean_POD20 = mean(error_POD20,1);
    mean_POD50 = mean(error_POD50,1);
    mean_POD1000 = mean(error_POD1000,1);
    
else 

    %% Calculate POD basis using 50 samples
    [~,~,Phi,~] = gen_samples(10);
    d_f_POD10=zeros(1,m);
    error_POD10 = zeros(N_rep,m);

end

%% Adjoint 
P = zeros(nele,m);
P_noise = zeros(nele,20);
for i=1:20
    z = zeros(m,1);
    if i<m+1
        z(i)=1;
        p = K'\(C'*z);
        P(:,i)=p;
    end
    z = randn(m,1);
    p = K'\(C'*z);
    P_noise (:,i)=p;
end
[psi,d,~]=svd(P);
W_adj = psi(:,1:m);
[W_noise,d2,~]=svd(P_noise);

d_f_LI=zeros(1,m);
d_f_OLR=zeros(1,m);
d_f_Sta=zeros(1,m);
d_f_POD_adj = zeros(1,m);
d_f_POD_noise = zeros(1,m);

for i=1:m

    %% LIS reduced operator
    [gamma_pos_LI,G_LI_cell{i},d_f_LI(i)]=solveReducedModel(V(:,1:i),W(:,1:i));
    [gamma_pos_Sta,G_Sta_cell{i},d_f_Sta(i)]=solveReducedModel(V_State(:,1:i),W(:,1:i));
    
    %% POD reduced operator
    if ~run_POD_analysis
        [gamma_pos_POD,G_POD_cell10{i},d_f_POD10(i)]=solveReducedModel(Phi(:,1:i),Phi(:,1:i));
        [~,G_POD_adj{i},d_f_POD_adj(i)]=solveReducedModel(Phi(:,1:i),W_adj(:,1:i));
        [~,G_POD_noise{i},d_f_POD_noise(i)]=solveReducedModel(Phi(:,1:i),W_noise(:,1:i));
    end

    %% Spantini Reduction
    [gamma_pos_OLR,G_OLR_cell{i},d_f_OLR(i)]=solveOLRA(V(:,1:i),W(:,1:i));
    
    %gamma_prior_red = Phi(:,1:i)'*gamma_prior_f*Phi(:,1:i);
    %G_red = C*Phi(:,1:i)*inv(Phi(:,1:i)'*K*Phi(:,1:i));
    %gamma_pos_red = gamma_prior_red-gamma_prior_red*G_red'*((G_red*gamma_prior_red*G_red'+gamma_obs)\G_red)*gamma_prior_red;
    %gamma_pos_PO = gamma_prior_f-Phi(:,1:i)*gamma_prior_red*G_red'*((G_red*gamma_prior_red*G_red'+gamma_obs)\G_red)*gamma_prior_red*Phi(:,1:i)';
    %gamma_pos_PO2 = gamma_prior_f-gamma_prior_f*Phi(:,1:i)*G_red'*((G_red*gamma_prior_red*G_red'+gamma_obs)\G_red)*Phi(:,1:i)'*gamma_prior_f;
    %d_f_POD2(i)=foerstnerDistance(gamma_pos_PO);
    
end

%% Posterior Mean Analysis

error_LI = zeros(N_rep,m);
error_OLR = zeros(N_rep,m);
error_Sta = zeros(N_rep,m);

error_adj = zeros(N_rep,m);
error_noise = zeros(N_rep,m);

mu_tilde = meanCalculation(G,zeros(m,1));

for j=1:N_rep
    ysam = C*state_sample(:,j)+sqrt(gamma_obs)*randn(m,1);
    % Full model mean
    mu_full = meanCalculation(G,ysam);
    mu_full_norm = norm(mu_full);

    for i = 1:m
        % LIS approximation
        mu_LIS = meanCalculation(G_LI_cell{i},ysam);
        mu_Sta = meanCalculation(G_Sta_cell{i},ysam);

        % POD approximation
        if ~run_POD_analysis
            mu_POD = meanCalculationPOD(G_POD_cell10{i},ysam,Phi(:,1:i));
            error_POD10(j,i) = norm(mu_full-mu_POD)/mu_full_norm;
            mu_POD = meanCalculation(G_POD_adj{i},ysam);
            error_adj(j,i) = norm(mu_full-mu_POD)/mu_full_norm;
            mu_POD = meanCalculation(G_POD_noise{i},ysam);
            error_noise(j,i) = norm(mu_full-mu_POD)/mu_full_norm;
        end

        % Spantini approximation
        mu_Sp_3 = meanCalculation(G_OLR_cell{i},ysam);

        error_LI(j,i) = norm(mu_full-mu_LIS)/mu_full_norm;
        error_OLR(j,i) = norm(mu_full-mu_Sp_3)/mu_full_norm;
        error_Sta(j,i) = norm(mu_full-mu_Sta)/mu_full_norm;
        
    end
    
end

mean_LI = mean(error_LI,1);
if ~run_POD_analysis
    mean_POD10 = mean(error_POD10,1);
end
mean_OLR = mean(error_OLR,1);

mean_Sta = mean(error_Sta,1);

mean_POD_adj = mean(error_adj,1);
mean_POD_noise = mean(error_noise,1);

%% PLOTS

width = 14.5;
height= 8;

if run_POD_analysis
    
    figure
    t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    ax1 = nexttile;
    semilogy(mean_POD10,'LineWidth',2)
    set(gca,'FontSize',20)
    box off
    hold on
    semilogy(mean_POD20,'LineWidth',2,'LineStyle','--')
    semilogy(mean_POD50,'LineWidth',2,'LineStyle',':')
    semilogy(mean_POD1000,'LineWidth',2,'LineStyle','-.')
    legend('10','20','50','1000','Location','southwest')
    legend boxoff
    title('Relative posterior mean error','Interpreter','latex','FontSize',28)
    axis([1 10 1e-5 .1])
    yticks([1e-5 1e-4 1e-3 1e-2 1e-1])
    xlabel('Approximation rank $r$','Interpreter','latex')
    ax=gca;
    
    % Second plot
    ax = nexttile;
    semilogy(sqrt(d_f_POD10),'LineWidth',2)
    set(gca,'FontSize',20)
    box off
    hold on
    semilogy(sqrt(d_f_POD20),'LineWidth',2,'LineStyle','--')
    semilogy(sqrt(d_f_POD50),'LineWidth',2,'LineStyle',':')
    semilogy(sqrt(d_f_POD1000),'LineWidth',2,'LineStyle','-.')
    title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',28)
    %ylabel('$d_F$','Interpreter','Latex',FontSize=14)
    xlabel('Approximation rank $r$','Interpreter','latex')
    legend('10','20','50','1000','Location','southwest')
    legend boxoff
    axis([1 10 1e-5 .1])
    yticks([1e-5 1e-4 1e-3 1e-2 1e-1])
    
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [.5 .5 width height]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [width height]);
    set(gcf, 'PaperPosition', [0 0 width height]);
    exportgraphics(gcf, 'PODBar.pdf', 'ContentType', 'vector');
end

figure
t = tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

% First plot
ax1 = nexttile;
imagesc(x_dofs,x_dofs,gamma_prior_f)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ylabel('z','Interpreter','latex')
title("Prior covariance",'Interpreter','latex','FontSize',22)
n = 256;
cmap = [ ...
    linspace(1,0,n)', ...   % R: 1 → 0
    linspace(1,0,n)', ...   % G: 1 → 0
    ones(n,1) ];           % B: stays 1

colormap(cmap)

ax2 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
%title('Analytical posterior covariance $\mathbf{\Gamma}_{\mathrm{pos}}$','Interpreter','latex','FontSize',18)
title('Posterior covariance','Interpreter','latex','FontSize',22)
colormap(cmap)

ax3 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos_LI)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
colorbar
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\scriptscriptstyle \mathrm{LIS}}$', ...
      'Interpreter', 'latex', 'FontSize', 22)
cb = colorbar(ax3, 'Location', 'eastoutside');
colormap(cmap)

nexttile;
axis off;

ax5 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos_OLR)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ylabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
colormap(cmap)
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\scriptscriptstyle \mathrm{OLR}}$','Interpreter','latex','FontSize',22)

ax6 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos_POD)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
colorbar
colormap(cmap)
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\scriptscriptstyle \mathrm{POD}}$','Interpreter','latex','FontSize',22)

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

%exportgraphics(gcf, 'priorBar.pdf', 'ContentType', 'vector');

height = 6;
alpha = 0.25;
LI_color = (1-alpha)*[0.4660 0.6740 0.1880]+alpha*[1 1 1];
OLR_color = (1-alpha)*[0.8500 0.3250 0.0980]+alpha*[1 1 1];
alpha=0.0;
POD_color = (1-alpha)*[0.3010 0.7450 0.9330]+alpha*[1 1 1];

figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
semilogy(mean_LI,'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',20)
box off
hold on
semilogy(mean_POD10,'--','Color',POD_color,'LineWidth',2)
semilogy(mean_OLR,'o','Color',OLR_color,'LineWidth',2)
semilogy(mean_POD_adj,"LineWidth",2)
semilogy(mean_POD_noise,"r","LineWidth",2)
%semilogy(mean_Sta,'LineWidth',2)

legend('LIS','POD','OLR',"adjoint","noise",'Location','southwest')
legend boxoff
title('Relative posterior mean error','Interpreter','latex','FontSize',28)
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])
xlabel('Approximation rank $r$','Interpreter','latex')

% Second plot
ax = nexttile;
semilogy(sqrt(d_f_LI),'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',20)
box off
hold on
semilogy(sqrt(d_f_POD10),'--','Color',POD_color,'LineWidth',2)
semilogy(sqrt(d_f_OLR),'o','Color',OLR_color,'LineWidth',2)
semilogy(sqrt(d_f_POD_adj),"LineWidth",2)
semilogy(sqrt(d_f_POD_noise),"r","LineWidth",2)
%semilogy(sqrt(d_f_Sta),'LineWidth',2)
title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',28)
%ylabel('F$\ddot{o}$rstner distance','Interpreter','latex')
xlabel('Approximation rank $r$','Interpreter','latex')
legend('LIS','POD','OLR',"Adjoint", "Adj noise",'Location','southwest')
legend boxoff
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'posBar.pdf', 'ContentType', 'vector');