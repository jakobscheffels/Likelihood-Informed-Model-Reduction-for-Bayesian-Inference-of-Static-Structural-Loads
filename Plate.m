clear all
addpath src\

rng(42)
load plate_no_hole.mat

nodes = size(coordinatesFEM,1);
num_ele = size(elementsFEM,1);
ndof = 2*nodes;

Y_bottom = min(coordinatesFEM(:,2));
tol = 1e-8;
plate_width = max(coordinatesFEM(:,1))-min(coordinatesFEM(:,1));

bottom_nodes = find(abs(coordinatesFEM(:,2) - Y_bottom) < tol);
dofs_bottom = zeros(2*size(bottom_nodes,1),1);
for i=1:size(bottom_nodes,1)
    dofs_bottom(2*i-1) = 2.*bottom_nodes(i)-1;
    dofs_bottom(2*i) = 2.*bottom_nodes(i);
end
coords_bottom = coordinatesFEM(bottom_nodes,:);

Y_top = max(coordinatesFEM(:,2));
top_nodes = find(abs(coordinatesFEM(:,2)-Y_top)<tol);
dofs_top = zeros(2*size(top_nodes,1),1);
dofs_topY = zeros(size(top_nodes,1),1);
for i=1:size(top_nodes,1)
    dofs_top(2*i-1)=2.*top_nodes(i)-1;
    dofs_top(2*i)=2.*top_nodes(i);
    dofs_topY(i)=2.*top_nodes(i);
end
coords_top = coordinatesFEM(top_nodes,:);
plate_height = coords_top(1,2)-coords_bottom(1,2);


%% Loading 
l_c = 0.3*plate_width;
mu=1e6;
delta=0.3;
mu_q = zeros(length(top_nodes),1);

nSegments = length(top_nodes)-1;
len = zeros(nSegments,1);
x_mid = zeros(nSegments,1);
for i=1:nSegments
    x_left = coordinatesFEM(top_nodes(i),1);
    x_right = coordinatesFEM(top_nodes(i+1),1);
    x_mid(i) = (x_right+x_left)/2;
    len(i)=x_right-x_left;
end
sigma = delta*mu;
B = zeros(nSegments,nSegments);
for i=1:nSegments
    for j=1:nSegments
        B(i,j)=sigma^2*exp(-abs(x_mid(i)-x_mid(j))/l_c);
    end
    mu_q(i:i+1)=mu_q(i:i+1)+len(i)/2.*[mu;mu];
end
L = chol(B,"lower");
z=randn(nSegments,1);
q_realization = mu+L*z;

mu_F = zeros(ndof,1);
F = zeros(ndof,1);

q_POD = mu+L*randn(nSegments,10);
F_POD = zeros(ndof,10);

N_rep = 200;
q_mean = mu+L*randn(nSegments,N_rep);
F_mean = zeros(ndof,N_rep);


for i=1:nSegments
    node_left = top_nodes(i);
    node_right = top_nodes(i+1);
    L_seg = coordinatesFEM(node_right,1)-coordinatesFEM(node_left,1);
    f_local = q_realization(i)*L_seg/2 * [1;1];

    F(dofs_top(2*i))=F(dofs_top(2*i))+f_local(1);
    F(dofs_top(2*(i+1)))=F(dofs_top(2*(i+1)))+f_local(2);

    f_local = q_POD(i,:).*L_seg/2 .* [1;1];
    F_POD(dofs_top(2*i),:)=F_POD(dofs_top(2*i),:)+f_local(1,:);
    F_POD(dofs_top(2*(i+1)),:)=F_POD(dofs_top(2*(i+1)),:)+f_local(2,:);

    f_local = q_mean(i,:).*L_seg/2.*[1;1];
    F_mean(dofs_top(2*i),:)=F_mean(dofs_top(2*i),:)+f_local(1,:);
    F_mean(dofs_top(2*(i+1)),:)=F_mean(dofs_top(2*(i+1)),:)+f_local(2,:);
end
mu_F(dofs_topY)=mu_q;
% for plotting
F_full = F;


%% IP Parameters
num_col = 2;
num_rows = 3;
gamma_obs = 1e-4;

m = num_col*num_rows;
m=6;

hx = plate_width/(num_col+1);
hy = plate_height/(num_rows+1);
m_pos = zeros(m,2);
m_dofs = zeros(m,1);
for i=1:num_col
    for j=1:num_rows
        posX = coords_bottom(1,1)+i*hx;
        posY = coords_bottom(1,2)+j*hy;
        [idx,D]=knnsearch(coordinatesFEM,[posX posY]);
        m_pos((i-1)*num_rows+j,:)=coordinatesFEM(idx,:);
        m_dofs((i-1)*num_rows+j)=idx;
    end
end

% at top 
m_pos = randi([1,17],1,m);
m_pos = unique(m_pos);

while ~(size(unique(m_pos),2)==m)
    m_pos=[m_pos, randi([1,17],1,m-size(m_pos,2))];
    m_pos = unique(m_pos);
end

m_pos=sort(m_pos);
m_dofs = dofs_topY(m_pos);
m_dofs = sort(m_dofs);

gamma_obs = gamma_obs.^2*eye(m);
S_obs = sqrt(gamma_obs);

gamma_prior_q=zeros(nSegments,nSegments);
x_coords_top = coords_top(:,1);
for i = 1:nSegments
    for j=1:nSegments
        gamma_prior_q(i,j)=sigma^2*exp(-abs((x_coords_top(i)+x_coords_top(i+1))/2-((x_coords_top(j)+x_coords_top(j+1))/2))/l_c);
    end
end
S_q = chol(gamma_prior_q,"lower");
S_pr = zeros(ndof,size(S_q,1));
S_pr (dofs_topY(1:end-1),:)=L_seg/2.*S_q;
S_pr(dofs_topY(2:end),:) =S_pr(dofs_topY(2:end),:)+L_seg/2.*S_q;
gamma_prior = zeros(ndof,ndof);
gamma_prior(dofs_topY(1:end-1),dofs_topY(1:end-1))=L_seg/2.*gamma_prior_q;
gamma_prior(dofs_topY(2:end),dofs_topY(2:end))=gamma_prior(dofs_topY(2:end),dofs_topY(2:end))+L_seg/2.*gamma_prior_q;
gamma_prior_full = gamma_prior;
gamma_prior(dofs_bottom,:)=[];
gamma_prior(:,dofs_bottom)=[];
S_pr (dofs_bottom,:)=[];

%m_pos = randi(nodes,m,1);
%m_pos = sort(m_pos);
C = zeros(m,ndof);
%C(1:m,2.*sort(m_dofs))=eye(m);
C (:,m_dofs)=eye(m);
C (:,dofs_bottom)=[];

%% --- Load mesh ---
% Mesh file should contain:
% nodes      -> N x 2 array of [x y] coordinates
% elements   -> M x 4 array of node indices for quads
%load('plate_hole.mat');  

Nnodes = size(coordinatesFEM,1);
Nelements = size(elementsFEM,1);
ndof = 2*Nnodes;

%% Elemental degree of freedom table
elementDofs = zeros(Nelements,8);
for i = 1:Nelements
    idx = elementsFEM(i,:);
    for j=1:4
        elementDofs(i,2*j-1:2*j)=[2*idx(j)-1 2*idx(j)];
    end
end

%% --- Material properties ---
E  = 210e9;  % Young's modulus [Pa]
nu = 0.3;    % Poisson's ratio

%% --- Gauss points for 2x2 integration ---
K = zeros(ndof,ndof);
for i=1:Nelements
    Ke = elementalStiffnessPlate(coordinatesFEM(elementsFEM(i,:)',:),E,nu);
    idx = elementDofs(i,:);
    K(idx,idx)=K(idx,idx)+Ke;
end

%% Apply Boundary Condition
K(dofs_bottom,:)=[];
K(:,dofs_bottom)=[];
F(dofs_bottom)=[];
mu_F(dofs_bottom)=[];

F_POD(dofs_bottom,:)=[];
F_mean(dofs_bottom,:)=[];

G = K\C';
G = G';

gamma_pos = gamma_prior-gamma_prior*G'*((G*gamma_prior*G'+gamma_obs)\G)*gamma_prior';

%% --- Solve ---
U = K\F;
Y = C*U+S_obs*randn(m,1);
G=C/K;

Y_mean = G*F_mean+S_obs*randn(m,N_rep);

mu_pos = mu_F+ gamma_pos*G'*(gamma_obs\(Y-G*mu_F));

U_full = zeros(ndof,1);
freeDofs = (1:1:ndof)';
freeDofs(dofs_bottom)=[];
[~,loc]=ismember(dofs_topY,freeDofs);

U_full(freeDofs)=U;

top_v_free = intersect(freeDofs,dofs_topY);

%% Approximations 

[psi,~,~]=svd(K\F_POD);

[phi,delta,nu]=svd((S_obs\G)*S_pr);
%delta=diag(delta);
%V = K*C'*((C*C')\S_obs)*phi(:,1:m)*delta(1:m,1:m);
V = S_pr*nu(:,1:m);
V_St = C'*((C*C')\S_obs)*phi(:,1:m);
W = G'*((S_obs\phi(:,1:m))/delta(1:m,1:m));
%V = gamma_prior*W;
delta = diag(delta);

% OLR
d_f_Sp = zeros(1,m);
p_OLR = zeros(1,m);

for r=1:m
    G_r = G*V(:,1:r)*W(:,1:r)';
    G_OLR_cell{r}=G_r;
    gamma_pos_Sp = gamma_prior-gamma_prior*G_r'*((G_r*gamma_prior*G_r'+gamma_obs)\G_r)*gamma_prior';
    d_f_Sp(r)=foerstnerDistance(gamma_pos_Sp,gamma_pos,W);
    p_OLR(r)=norm(gamma_pos-gamma_pos_Sp,2);
end

% LIS
d_f_LIS = zeros(1,m);
p_LIS = zeros(1,m);
pHessLIS = eps.*ones(1,m);
bound_LIS = eps.*ones(1,m);
C_LIS = zeros(1,m);
for r=1:m
    K_hat = W(:,1:r)'*K*V(:,1:r);
    G_hat = C*V(:,1:r)*(K_hat\W(:,1:r)');
    G_LIS_cell{r}=G_hat;
    gamma_pos_LIS = gamma_prior-gamma_prior*G_hat'*((G_hat*gamma_prior*G_hat'+gamma_obs)\G_hat)*gamma_prior';
    d_f_LIS(r)=foerstnerDistance(gamma_pos_LIS,gamma_pos,W);
    p_LIS(r)=norm(gamma_pos_LIS-gamma_pos,2);
    pHessLIS(r)=norm((S_obs\(G-G_hat))*S_pr,2);
    if r<m
        R = phi(:,1:r)'*(S_obs\C)*V(:,1:r);
        bound_LIS(r)=delta(r+1)+norm((phi(:,1:r)-(S_obs\C)*(V(:,1:r)/R))*diag(delta(1:r)),2);
    end
    C_LIS(r)=norm(S_pr,2)^2*(norm(((eye(m)+(S_obs\G)*gamma_prior*(G'/S_obs))\(S_obs\G)*S_pr),2) ...
        +norm((S_obs\G_hat)*S_pr,2)*norm((S_obs\G)*S_pr,2)*(norm((S_obs\G)*S_pr,2)+norm((S_obs\G_hat)*S_pr,2)) ...
        +norm(((eye(m)+(S_obs\G_hat)*gamma_prior*(G_hat'/S_obs))\(S_obs\G_hat)*S_pr),2));
end

d_f_Sta = zeros(1,m);
for r=1:m
    K_hat = W(:,1:r)'*K*V_St(:,1:r);
    G_hat = C*V_St(:,1:r)*(K_hat\W(:,1:r)');
    G_Sta_cell{r}=G_hat;
    gamma_pos_LIS = gamma_prior-gamma_prior*G_hat'*((G_hat*gamma_prior*G_hat'+gamma_obs)\G_hat)*gamma_prior';
    d_f_Sta(r)=foerstnerDistance(gamma_pos_LIS,gamma_pos,W);
end

% POD 
d_f_POD = zeros(1,m);
for r=1:m
    K_hat = psi(:,1:r)'*K*psi(:,1:r);
    G_hat = C*psi(:,1:r)*(K_hat\psi(:,1:r)');
    G_POD_cell{r}=G_hat;
    gamma_pos_POD = gamma_prior-gamma_prior*G_hat'*((G_hat*gamma_prior*G_hat'+gamma_obs)\G_hat)*gamma_prior';
    d_f_POD(r)=foerstnerDistance(gamma_pos_POD,gamma_pos,W);
end

%% Posterior mean
error_LI=zeros(N_rep,m);
error_POD=zeros(N_rep,m);
error_OLR=zeros(N_rep,m);
error_Sta = zeros(N_rep,m);

for j=1:N_rep
    mu_full = mu_F+ gamma_pos*G'*(gamma_obs\(Y_mean(:,j)-G*mu_F));
    mu_full_norm = norm(mu_full);
    for r=1:m
        G_hat = G_LIS_cell{r};
        mu_LI = mu_F+gamma_prior*G_hat'*((G*gamma_prior*G'+gamma_obs)\(Y_mean(:,j)-G_hat*mu_F));
        G_hat = G_Sta_cell{r};
        mu_Sta = mu_F+gamma_prior*G_hat'*((G*gamma_prior*G'+gamma_obs)\(Y_mean(:,j)-G_hat*mu_F));
        G_hat = G_POD_cell{r};
        mu_POD = mu_F+gamma_prior*G_hat'*((G*gamma_prior*G'+gamma_obs)\(Y_mean(:,j)-G_hat*mu_F));
        G_hat = G_OLR_cell{r};
        mu_OLR = mu_F+gamma_prior*G_hat'*((G*gamma_prior*G'+gamma_obs)\(Y_mean(:,j)-G_hat*mu_F));

        error_LI(j,r) = norm(mu_full-mu_LI)/mu_full_norm;
        error_POD(j,r) = norm(mu_full-mu_POD)/mu_full_norm;
        error_OLR(j,r) = norm(mu_full-mu_OLR)/mu_full_norm;
        error_Sta(j,r) = norm(mu_full-mu_Sta)/mu_full_norm;
    end
end

mean_LI = mean(error_LI,1);
mean_POD = mean(error_POD,1);
mean_OLR = mean(error_OLR,1);
mean_Sta = mean(error_Sta,1);

%% --- Plot deformed mesh ---
scale=1000;

%% Posterior plots
height= 8;
width = 14.5;
alpha = 0.25;
LI_color = (1-alpha)*[0.4660 0.6740 0.1880]+alpha*[1 1 1];
OLR_color = (1-alpha)*[0.8500 0.3250 0.0980]+alpha*[1 1 1];
alpha=0.00;
POD_color = (1-alpha)*[0.3010 0.7450 0.9330]+alpha*[1 1 1];

figure
%tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
%nexttile;
semilogy(mean_LI,'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',24)
box off
hold on
semilogy(mean_POD,'--','Color',POD_color,'LineWidth',2)
semilogy(mean_OLR,'o','Color',OLR_color,'LineWidth',2)
%semilogy(mean_Sta,'LineWidth',2)

legend('LIS','POD','OLR','State','Location','southwest')
legend boxoff
title('Relative posterior mean error','Interpreter','latex','FontSize',36)
axis([1 m 1e-15 100])
yticks([10^(-15) 10^(-10) 10^(-5) 10^0])
xlabel('Approximation rank $r$','Interpreter','latex')

% Second plot
%ax = nexttile;
figure
semilogy(sqrt(d_f_LIS),'Color',LI_color,'LineWidth',2)
set(gca,'FontSize',24)
box off
hold on
semilogy(sqrt(d_f_POD),'--','Color',POD_color,'LineWidth',2)
semilogy(sqrt(d_f_Sp),'o','Color',OLR_color,'LineWidth',2)
%semilogy(sqrt(d_f_Sta),'LineWidth',2)
%title('F$\ddot{o}$rstner posterior covariance error','Interpreter','latex','FontSize',36)
title('Posterior covariance error','Interpreter','latex','FontSize',36)
%ylabel('F$\ddot{o}$rstner distance','Interpreter','latex')
xlabel('Approximation rank $r$','Interpreter','latex')
legend('LIS','POD','OLR','State','Location','southwest')
legend boxoff
axis([1 m 1e-15 100])
yticks([10^(-15) 10^(-10) 10^(-5) 1])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [.5 .5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

%% BOUNDS

figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
nexttile
semilogy(bound_LIS,"LineWidth",2)
set(gca,"FontSize",20)
hold on
semilogy(pHessLIS,"k--","LineWidth",2)
legend("Bound","Actual","Interpreter","Latex","Location","southwest")
title("Preconditioned Hessian bound $p=2$","Interpreter","Latex","FontSize",20)

nexttile;
semilogy(C_LIS.*bound_LIS,"LineWidth",2)
set(gca,"FontSize",20)
hold on
semilogy(p_LIS,"k--","LineWidth",2)
legend("Bound","Actual","Interpreter","Latex","Location","southwest")
title("Posterior covariance error for $p=2$","Interpreter","latex","FontSize",20)


%% Structure
% Plot mesh and force vector

figure
hold on

for e = 1:size(elementsFEM,1)
%for e = 1:1
    idx = elementsFEM(e,:);
    
    x = coordinatesFEM(idx,1);
    y = coordinatesFEM(idx,2);
    
    % close the element polygon
    x = [x; x(1)];
    y = [y; y(1)];
    
    plot(x,y,'k')
end

for e = 1:size(elementsFEM,1)
    idx = elementsFEM(e,:);
    dofs = elementDofs(e,:);
    xdofs = dofs(1:2:end);
    ydofs = dofs(2:2:end);
    
    x = coordinatesFEM(idx,1)+scale.*U_full(xdofs);
    y = coordinatesFEM(idx,2)+scale.*U_full(ydofs);
    
    % close the element polygon
    x = [x; x(1)];
    y = [y; y(1)];
    
    plot(x,y,'r')
end



figure
hold on
rectangle("Position",[coords_bottom(1,1),coords_bottom(1,2),coords_bottom(end,1)-coords_bottom(1,1),coords_top(1,2)-coords_bottom(1,2)] ...
    ,"FaceColor",[.75 .75 .75]);
set(gca,"FontSize",24)
for e = 1:size(elementsFEM,1)
%for e = 1:1
    idx = elementsFEM(e,:);
    
    x = coordinatesFEM(idx,1);
    y = coordinatesFEM(idx,2);
    
    % close the element polygon
    x = [x; x(1)];
    y = [y; y(1)];
    
    plot(x,y,'k')
end
for i = 1:size(top_nodes,1)
    quiver(coords_top(i,1),coords_top(i,2),0,F_full(dofs_topY(i))./mu*5,"r","LineWidth",1.5)
end

%for i=1:m
%    plot(m_pos(i,1),m_pos(i,2),'bo',"MarkerSize",6,"LineWidth",1.5)
%end
for i=1:m
    plot(coords_top(m_pos(i),1),coords_top(m_pos(i),2),"bo","MarkerSize",6,"LineWidth",1.5)
end
xlabel("Width","Interpreter","latex")
ylabel("Height","Interpreter","latex")


figure
%tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
%nexttile
imagesc(gamma_prior(loc,loc))
set(gca,"FontSize",24)
xlabel("Width","Interpreter","latex")
title("Prior covariance","Interpreter","Latex","FontSize",38)

%ax1=nexttile;
figure
imagesc(gamma_pos(loc,loc))
set(gca,"FontSize",24)
xlabel("Width","Interpreter","latex")
title("Posterior covariance","Interpreter","Latex","FontSize",38)
colorbar
ax1.CLim=[0 max(max(gamma_prior))];

figure
plot(mu_F(loc))
set(gca,"FontSize",20)
hold on
plot(mu_pos(loc))


function df = foerstnerDistance(gamma1,gamma_pos,W)
    %load LIP_Setup.mat gamma_pos
    %load LIS_Basis.mat W
  
    % Calculate cholesky factor of projected covariances
    L = chol(W'*gamma1*W,'lower');
    R = chol(W'*gamma_pos*W,'lower');
    %L = chol(gamma1,'lower');
    %R = chol(gamma_pos,'lower');
    % Solve for generalized eigenvalues
    df=gen_eigenvalue(L,R);
    df(abs(df)<eps)=[];
    % Evaluate squared Foerstner distance
    df = dot(log(df),log(df));
end

function delta = gen_eigenvalue(L,R)
    % Solve for generalized eigenvalues using SVD
    tmp = L'/R';
    [~,S,~]=svd(tmp);
    delta = diag(S).^2;
end