%clearvars

addpath src\
%% Parameters of the problem
rng(25);
run_POD_analysis = false; % set 'true' if want to run POD analysis

%% Setup of the problem
Parameters({'beam_bool'},{true});
load('Parameters.mat')

x_dofs=applyBoundaryCondition(x_dofs,BC_dofs,'Coordinate');

LIP_Setup({'gamma_obs'},{0.001^2});
load('LIP_Setup.mat')

%% Generalized eigenvectors / LIS basis
%[~,~,omega]=svd((chol(gamma_obs,'lower')\G)*S_pr);
%V = S_pr*omega;
[delta,V,W,V_State] = calculateLISBasis();

%[~,state_samples]=gen_samples(100);
%save LIS_Basis_Beam.mat V W K mu_f state_samples

%% Auxiliary IP
[d,V_tilde,W_tilde]=calculateLISBasisIP();

%% Mean samples
N_rep = 200;

% Generate sample
[~,state_sample]=gen_samples(N_rep);


%% Calculate POD basis using 50 samples
[~,~,Phi,~] = gen_samples(10);
d_f_POD10=zeros(1,m);
error_POD10 = zeros(N_rep,m);

d_f_LI=zeros(1,m);
d_f_OLR=zeros(1,m);
d_f_Sta=zeros(1,m);
d_f_IP = zeros(1,m);

for i=1:m

    %% LIS reduced operator
    [gamma_pos_LI,G_LI_cell{i},d_f_LI(i)]=solveReducedModel(V(:,1:i),W(:,1:i));
    [gamma_pos_Sta,G_Sta_cell{i},d_f_Sta(i)]=solveReducedModel(V_State(:,1:i),W(:,1:i));
    [gamma_pos_IP,G_IP_cell{i},d_f_IP(i)]=solveReducedModel(V_tilde(:,1:i),W_tilde(:,1:i));
    
    %% POD reduced operator
    if ~run_POD_analysis
        [gamma_pos_POD,G_POD_cell10{i},d_f_POD10(i)]=solveReducedModel(Phi(:,1:i),Phi(:,1:i));
    end

    %% Spantini Reduction
    [gamma_pos_OLR,G_OLR_cell{i},d_f_OLR(i)]=solveOLRA(V(:,1:i),W(:,1:i));
    
    gamma_prior_red = Phi(:,1:i)'*gamma_prior_f*Phi(:,1:i);
    G_red = C*Phi(:,1:i)*inv(Phi(:,1:i)'*K*Phi(:,1:i));
    gamma_pos_red = gamma_prior_red-gamma_prior_red*G_red'*((G_red*gamma_prior_red*G_red'+gamma_obs)\G_red)*gamma_prior_red;
    %gamma_pos_PO = gamma_prior_f-Phi(:,1:i)*gamma_prior_red*G_red'*((G_red*gamma_prior_red*G_red'+gamma_obs)\G_red)*gamma_prior_red*Phi(:,1:i)';
    %gamma_pos_PO2 = gamma_prior_f-gamma_prior_f*Phi(:,1:i)*G_red'*((G_red*gamma_prior_red*G_red'+gamma_obs)\G_red)*Phi(:,1:i)'*gamma_prior_f;
    %d_f_POD2(i)=foerstnerDistance(gamma_pos_PO);
    
end

%% Posterior Mean Analysis

error_LI = zeros(N_rep,m);
error_OLR = zeros(N_rep,m);
error_Sta = zeros(N_rep,m);
error_IP = zeros(N_rep,m);

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
        mu_IP = meanCalculation(G_IP_cell{i},ysam);

        % POD approximation
        if ~run_POD_analysis
            mu_POD = meanCalculation(G_POD_cell10{i},ysam);
            error_POD10(j,i) = norm(mu_full-mu_POD)/mu_full_norm;
        end

        % Spantini approximation
        mu_Sp_3 = meanCalculation(G_OLR_cell{i},ysam);

        error_LI(j,i) = norm(mu_full-mu_LIS)/mu_full_norm;
        error_OLR(j,i) = norm(mu_full-mu_Sp_3)/mu_full_norm;
        error_Sta(j,i) = norm(mu_full-mu_Sta)/mu_full_norm;
        error_IP(j,i) = norm(mu_full-mu_IP)/mu_full_norm;
    end
    
end

mean_LI = mean(error_LI,1);
if ~run_POD_analysis
    mean_POD10 = mean(error_POD10,1);
end
mean_OLR = mean(error_OLR,1);

mean_Sta = mean(error_Sta,1);
mean_IP = mean(error_IP,1);

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
semilogy(mean_Sta,'LineWidth',2)
semilogy(mean_IP,'r--','LineWidth',2)
legend('LIS','POD','OLR','SVD','IP','Location','southwest')
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
semilogy(sqrt(d_f_Sta),'LineWidth',2)
semilogy(sqrt(d_f_IP),'r--','LineWidth',2)
title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',28)
%ylabel('F$\ddot{o}$rstner distance','Interpreter','latex')
xlabel('Approximation rank $r$','Interpreter','latex')
legend('LIS','POD','OLR','SVD','IP','Location','southwest')
legend boxoff
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

figure
semilogy(Delta,'LineWidth',2)
hold on
%semilogy(Delta2,"LineWidth",2)
semilogy(delta,'LineWidth',2)
title("Singular values","Interpreter","latex","FontSize",22)
%legend("$\gamma=10^3$","$\gamma=10^7$","$\Delta$","Interpreter","latex")

%exportgraphics(gcf, 'posBar.pdf', 'ContentType', 'vector');