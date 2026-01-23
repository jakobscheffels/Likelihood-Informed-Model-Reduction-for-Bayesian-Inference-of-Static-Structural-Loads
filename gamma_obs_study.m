load Bar_obs.mat
mean_LI5 = mean_LI;
mean_POD5 = mean_POD10;
mean_OLR5 = mean_OLR;
d_f_LI5 = d_f_LI;
d_f_POD5 = d_f_POD10;
d_f_OLR5 = d_f_OLR;

load Bar_obs0001.mat
mean_LI3 = mean_LI;
mean_POD3 = mean_POD10;
mean_OLR3 = mean_OLR;
d_f_LI3 = d_f_LI;
d_f_POD3 = d_f_POD10;
d_f_OLR3 = d_f_OLR;

load Bar_obs01.mat
mean_LI1 = mean_LI;
mean_POD1 = mean_POD10;
mean_OLR1 = mean_OLR;
d_f_LI1 = d_f_LI;
d_f_POD1 = d_f_POD10;
d_f_OLR1 = d_f_OLR;

load Bar_obs001.mat
mean_LI2 = mean_LI;
mean_POD2 = mean_POD10;
mean_OLR2 = mean_OLR;
d_f_LI2 = d_f_LI;
d_f_POD2 = d_f_POD10;
d_f_OLR2= d_f_OLR;

load Bar_obs00001.mat
mean_LI4 = mean_LI;
mean_POD4 = mean_POD10;
mean_OLR4 = mean_OLR;
d_f_LI4 = d_f_LI;
d_f_POD4 = d_f_POD10;
d_f_OLR4 = d_f_OLR;

load Bar_obs7.mat
mean_LI7 = mean_LI;
mean_POD7 = mean_POD10;
mean_OLR7 = mean_OLR;
d_f_LI7 = d_f_LI;
d_f_POD7 = d_f_POD10;
d_f_OLR7 = d_f_OLR;

load Bar_obs6.mat
mean_LI6 = mean_LI;
mean_POD6 = mean_POD10;
mean_OLR6 = mean_OLR;
d_f_LI6 = d_f_LI;
d_f_POD6 = d_f_POD10;
d_f_OLR6 = d_f_OLR;

load Bar_obs1.mat
mean_LI = mean_LI;
mean_POD = mean_POD10;
mean_OLR = mean_OLR;

height = 6;
alpha = 0.5;
LI_color = (1-alpha)*[0.4660 0.6740 0.1880]+alpha*[1 1 1];
OLR_color = (1-alpha)*[0.8500 0.3250 0.0980]+alpha*[1 1 1];
alpha=0.25;
POD_color = (1-alpha)*[0.3010 0.7450 0.9330]+alpha*[1 1 1];

figure
t = tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
semilogy(mean_LI,'b','LineWidth',2)
set(gca,'FontSize',20)
box off
hold on
semilogy(mean_POD10,'b','LineWidth',2)
semilogy(mean_OLR,'bo','LineWidth',2)
semilogy(mean_LI3,'r','LineWidth',2)
semilogy(mean_POD3,'r','LineWidth',2)
semilogy(mean_OLR3,'ro','LineWidth',2)
semilogy(mean_LI5,'k','LineWidth',2)
semilogy(mean_POD5,'k','LineWidth',2)
semilogy(mean_OLR5,'ko','LineWidth',2)
legend('$10^{0}$','','','$10^{-3}$','','','$10^{-5}$','Location','southwest','Interpreter','latex')
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
semilogy(sqrt(d_f_POD10),'Color',POD_color,'LineWidth',2)
semilogy(sqrt(d_f_OLR),'o','Color',OLR_color,'LineWidth',2)
title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',28)
%ylabel('F$\ddot{o}$rstner distance','Interpreter','latex')
xlabel('Approximation rank $r$','Interpreter','latex')
legend('LIS','POD','OLR','Location','southwest')
legend boxoff
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

figure
t = tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = nexttile;
semilogy(mean_OLR)
set(gca,'FontSize',20)
hold on
semilogy(mean_OLR1)
semilogy(mean_OLR2)
semilogy(mean_OLR3)
semilogy(mean_OLR4)
semilogy(mean_OLR5)
semilogy(mean_OLR6)
semilogy(mean_OLR7,'--')
title("OLR",'Interpreter','latex','FontSize',22)
ylabel('Mean error')
xlabel('r','Interpreter','latex')
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

ax = nexttile;
semilogy(mean_LI)
set(gca,'FontSize',20)
hold on
semilogy(mean_LI1)
semilogy(mean_LI2)
semilogy(mean_LI3)
semilogy(mean_LI4)
semilogy(mean_LI5)
semilogy(mean_LI6)
semilogy(mean_LI7,'--')
title("LIS",'Interpreter','latex','FontSize',22)
xlabel('r','Interpreter','latex')
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

ax = nexttile;
semilogy(mean_POD)
set(gca,'FontSize',20)
hold on
semilogy(mean_POD1)
semilogy(mean_POD2)
semilogy(mean_POD3)
semilogy(mean_POD4)
semilogy(mean_POD5)
semilogy(mean_POD6)
semilogy(mean_POD7,'--')
xlabel('r','Interpreter','latex')
title("POD",'Interpreter','latex','FontSize',22)
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])
legend('$10^0$','$10^{-1}$','$10^{-2}$','$10^{-3}$','$10^{-4}$','$10^{-5}$','$10^{-6}$','$10^{-7}$','Interpreter','latex','Location','southwest','FontSize',14)
legend boxoff

ax = nexttile;
semilogy(sqrt(d_f_OLR))
set(gca,'FontSize',20)
hold on
semilogy(sqrt(d_f_OLR1))
semilogy(sqrt(d_f_OLR2))
semilogy(sqrt(d_f_OLR3))
semilogy(sqrt(d_f_OLR4))
semilogy(sqrt(d_f_OLR5))
semilogy(sqrt(d_f_OLR6))
semilogy(sqrt(d_f_OLR7),'--')
xlabel('r','Interpreter','latex')
ylabel('Cov error')
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

ax = nexttile;
semilogy(sqrt(d_f_LI))
set(gca,'FontSize',20)
hold on
semilogy(sqrt(d_f_LI1))
semilogy(sqrt(d_f_LI2))
semilogy(sqrt(d_f_LI3))
semilogy(sqrt(d_f_LI4))
semilogy(sqrt(d_f_LI5))
semilogy(sqrt(d_f_LI6))
semilogy(sqrt(d_f_LI7),'--')
xlabel('r','Interpreter','latex')
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

ax = nexttile;
semilogy(sqrt(d_f_POD10))
set(gca,'FontSize',20)
hold on
semilogy(sqrt(d_f_POD1))
semilogy(sqrt(d_f_POD2))
semilogy(sqrt(d_f_POD3))
semilogy(sqrt(d_f_POD4))
semilogy(sqrt(d_f_POD5))
semilogy(sqrt(d_f_POD6))
xlabel('r','Interpreter','latex')
semilogy(sqrt(d_f_POD7),'--')
axis([1 10 1e-18 1])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])


