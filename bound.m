load Bar.mat
S_obs = sqrt(gamma_obs);
[Omega,Delta,Nu]=svd((S_obs\G)*S_pr);
Delta=diag(Delta);
B = (S_obs\C)*S_pr*Nu;

sin_norm = zeros(1,10);
for r = 1:10
    B_r = B(:,1:r);
    Omega_r = Omega(:,1:r);
    [Q,R]=qr(B_r,0);
    M = Omega_r'*Q;
    [~,S,~]=svd(M);
    S=diag(S);
    sin_norm(r)=sqrt(1-S(end)^2);
end


error_bound = zeros(1,10);
%error_bound(end)=eps;
for i=1:9
    error_bound(i)=Delta(i+1)+Delta(1)*sin_norm(i);
end

%% Plots
width = 14.5;
height= 6;

figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
nexttile;
semilogy(Delta,'LineWidth',2)
set(gca,"FontSize",20)
axis([1 10 1e-8 1e2])
xlabel("r","Interpreter","latex")
title("Singular Values $\delta$ of LIS basis","FontSize",22, "Interpreter","latex")

nexttile;
semilogy(sin_norm,'LineWidth',2)
axis([1 10 1e-8 1e2])
set(gca,"FontSize",20)
xlabel("r","Interpreter","latex")
title("$\Vert sin \Theta(U_r,V_r)\Vert$","FontSize",22, "Interpreter","latex")

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [.5 .5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
nexttile;
semilogy(error_bound,'LineWidth',2)
set(gca,"FontSize",20)
axis([1 10 1e-14 10])
xlabel("r","Interpreter","latex")
title("Bound $\delta_{r+1}+\delta_1\cdot sin \Theta(U_r,V_r)$","FontSize",22,"Interpreter","latex")

nexttile;
semilogy(sqrt(d_f_LI),'LineWidth',2)
set(gca,"FontSize",20)
axis([1 10 1e-14 10])
xlabel("r","Interpreter","latex")
title("Covariance error","Interpreter","latex","FontSize",22)

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [.5 .5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

