clearvars

addpath src\
%% Parameters of the problem
rng(25);
run_POD_analysis = false; % set 'true' if want to run POD analysis

%% Setup of the problem
Parameters();
load('Parameters.mat')

x_dofs=applyBoundaryCondition(x_dofs,BC_dofs,'Coordinate');

LIP_Setup();
load('LIP_Setup.mat')

%% Generalized eigenvectors / LIS basis
[~,V,W] = calculateLISBasis();

%% Mean samples
N = 200;

% Generate sample
[~,state_sample]=gen_samples(N);

if run_POD_analysis
    d_f_POD10 = zeros(1,nm);
    d_f_POD20 = zeros(1,nm);
    d_f_POD50 = zeros(1,nm);
    d_f_POD1000 = zeros(1,nm);

    %% Calculate POD basis
    
    [~,~,Phi10,~] = gen_samples(10);
    [~,~,Phi20,~] = gen_samples(20);
    [~,~,Phi50,~] = gen_samples(50);
    [~,~,Phi1000,~] = gen_samples(1000);

    for i=1:nm
        [gamma_pos_POD,G_POD_cell50{i},d_f_POD50(i)]=solveReducedModel(Phi50(:,1:i),Phi50(:,1:i));

        [~,G_POD_cell10{i},d_f_POD10(i)]=solveReducedModel(Phi10(:,1:i),Phi10(:,1:i));
        [~,G_POD_cell20{i},d_f_POD20(i)]=solveReducedModel(Phi20(:,1:i),Phi20(:,1:i));
        [~,G_POD_cell1000{i},d_f_POD1000(i)]=solveReducedModel(Phi1000(:,1:i),Phi1000(:,1:i));
  
    end

    error_POD10 = zeros(N,nm);
    error_POD20 = zeros(N,nm);
    error_POD50 = zeros(N,nm);
    error_POD1000 = zeros(N,nm);

    for j=1:N
        ysam = B*state_sample(:,j)+sqrt(gamma_obs)*randn(nm,1);
        % Full model mean
        mu_full = meanCalculation(G,ysam);
        mu_full_norm = norm(mu_full);
    
        for i = 1:nm
    
            % POD approximation
            mu_POD50 = meanCalculation(G_POD_cell50{i},ysam);
            mu_POD10 = meanCalculation(G_POD_cell10{i},ysam);
            mu_POD20 = meanCalculation(G_POD_cell20{i},ysam);
            mu_POD1000 = meanCalculation(G_POD_cell1000{i},ysam);
  
  
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
    [~,~,Phi,~] = gen_samples(50);
    d_f_POD50=zeros(1,nm);
    error_POD50 = zeros(N,nm);

end

d_f_LI=zeros(1,nm);
d_f_OLR=zeros(1,nm);

for i=1:nm

    %% LIS reduced operator
    [gamma_pos_LI,G_LI_cell{i},d_f_LI(i)]=solveReducedModel(V(:,1:i),W(:,1:i));
    
    %% POD reduced operator
    if ~run_POD_analysis
        [gamma_pos_POD,G_POD_cell50{i},d_f_POD50(i)]=solveReducedModel(Phi(:,1:i),Phi(:,1:i));
    end

    %% Spantini Reduction
    [gamma_pos_OLR,G_OLR_cell{i},d_f_OLR(i)]=solveOLRA(V(:,1:i),W(:,1:i));
    
end

%% Posterior Mean Analysis

error_LI = zeros(N,nm);
error_OLR = zeros(N,nm);

mu_tilde = meanCalculation(G,zeros(nm,1));

for j=1:N
    ysam = B*state_sample(:,j)+sqrt(gamma_obs)*randn(nm,1);
    % Full model mean
    mu_full = meanCalculation(G,ysam);
    mu_full_norm = norm(mu_full);

    for i = 1:nm
        % LIS approximation
        mu_LIS = meanCalculation(G_LI_cell{i},ysam);

        % POD approximation
        if ~run_POD_analysis
            mu_POD = meanCalculation(G_POD_cell50{i},ysam);
            error_POD50(j,i) = norm(mu_full-mu_POD)/mu_full_norm;
        end

        % Spantini approximation
        mu_Sp_3 = meanCalculation(G_OLR_cell{i},ysam);

        error_LI(j,i) = norm(mu_full-mu_LIS)/mu_full_norm;
        error_OLR(j,i) = norm(mu_full-mu_Sp_3)/mu_full_norm;
    end
    
end

mean_LI = mean(error_LI,1);
if ~run_POD_analysis
    mean_POD50 = mean(error_POD50,1);
end
mean_OLR = mean(error_OLR,1);


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
    semilogy(mean_POD20,'LineWidth',2)
    semilogy(mean_POD50,'LineWidth',2)
    semilogy(mean_POD1000,'LineWidth',2)
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
    semilogy(sqrt(d_f_POD20),'LineWidth',2)
    semilogy(sqrt(d_f_POD50),'LineWidth',2)
    semilogy(sqrt(d_f_POD1000),'LineWidth',2)
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
title("Prior covariance",'Interpreter','latex','FontSize',28)

ax2 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
%title('Analytical posterior covariance $\mathbf{\Gamma}_{\mathrm{pos}}$','Interpreter','latex','FontSize',18)
title('Posterior covariance','Interpreter','latex','FontSize',28)

ax3 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos_LI)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
colorbar
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\scriptscriptstyle \mathrm{LI}}$', ...
      'Interpreter', 'latex', 'FontSize', 28)
cb = colorbar(ax3, 'Location', 'eastoutside');

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

title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\scriptscriptstyle \mathrm{OLR}}$','Interpreter','latex','FontSize',28)

ax6 = nexttile;
imagesc(x_dofs,x_dofs,gamma_pos_POD)
axis equal tight;
set(gca,'FontSize',20)
xlabel('z','Interpreter','latex')
ax = gca;
ax.CLim=[0 max(max(gamma_prior_f))];
yticks([])
colorbar
title('Approximation $\mathbf{\Gamma}_{\mathrm{pos}}^{\scriptscriptstyle \mathrm{POD}}$','Interpreter','latex','FontSize',28)

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

height = 6;
alpha = 0.5;
LI_color = (1-alpha)*[0.4660 0.6740 0.1880]+alpha*[1 1 1];
OLR_color = (1-alpha)*[0.8500 0.3250 0.0980]+alpha*[1 1 1];
alpha=0.25;
POD_color = (1-alpha)*[0.3010 0.7450 0.9330]+alpha*[1 1 1];

figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
semilogy(mean_LI,'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',20)
box off
hold on
semilogy(mean_POD50,'Color',POD_color,'LineWidth',2)
semilogy(mean_OLR,'o','Color',OLR_color,'LineWidth',2)

legend('LI','POD','OLR','Location','southwest')
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
semilogy(sqrt(d_f_POD50),'Color',POD_color,'LineWidth',2)
semilogy(sqrt(d_f_OLR),'o','Color',OLR_color,'LineWidth',2)
title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',28)
%ylabel('F$\ddot{o}$rstner distance','Interpreter','latex')
xlabel('Approximation rank $r$','Interpreter','latex')
legend('LI','POD','OLR','Location','southwest')
legend boxoff
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);