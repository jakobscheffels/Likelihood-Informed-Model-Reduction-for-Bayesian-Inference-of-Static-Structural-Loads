clear all

load plate_hole.mat

nodes = size(coordinatesFEM,1);
num_ele = size(elementsFEM,1);
ndof = 2*nodes;

Y_bottom = min(coordinatesFEM(:,2));
tol = 1e-8;

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
for i=1:size(top_nodes,1)
    dofs_top(2*i-1)=2.*top_nodes(i)-1;
    dofs_top(2*i)=2.*top_nodes(i);
end
coords_top = coordinatesFEM(top_nodes,:);

mu=1e6;
delta=0.03;
q=mu;

F = zeros(ndof,1);
F(dofs_top(2:2:end))=q.*(coords_top(2,1)-coords_top(1,1));
F(dofs_top(2))=F(dofs_top(2))/2;
F(top_nodes(end))=F(dofs_top(end))/2;

F_full = F;
%F(bottom_nodes)=[];

%E = 210e9;
%nu = 0.3;
%D= E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];


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

for e = 1:Nelements
    nodes_e = elementsFEM(e,:);      % 1x4 array of node indices
    coords_e = coordinatesFEM(nodes_e,:);  % 4x2 coordinates
    psi = 0; eta = 0;             % center of element
    J = Jacobian(coords_e(:,1), coords_e(:,2), eta, psi);
    if det(J) <= 0
        fprintf('Flipped element %d: det(J) = %f\n', e, det(J));
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

K(dofs_bottom,:)=[];
K(:,dofs_bottom)=[];
F(dofs_bottom)=[];

%% --- Solve ---
U = K\F;

U_full = zeros(ndof,1);
freeDofs = (1:1:ndof)';
freeDofs(dofs_bottom)=[];
U_full(freeDofs)=U;

%% --- Plot deformed mesh ---
scale=1000;

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
    %{
    list = [];
    listx = [];
    listy = [];
    for i=1:4
        if ismember(idx(i),bottom_nodes)
            list = [list i];
        end
        if ismember(xdofs(i),dofs_bottom)
            listx = [listx i];
        end
        if ismember(ydofs(i),dofs_bottom)
            listy = [listy i];
        end
    end
    idx(list)=[];
    xdofs(listx)=[];
    ydofs(listy)=[];
    %}
    x = coordinatesFEM(idx,1)+scale.*U_full(xdofs);
    y = coordinatesFEM(idx,2)+scale.*U_full(ydofs);
    
    % close the element polygon
    x = [x; x(1)];
    y = [y; y(1)];
    
    plot(x,y,'r')
end



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
for i = 1:size(top_nodes,1)
    quiver(coords_top(i,1),coords_top(i,2),0,F_full(dofs_top(2*i))./mu*5,"r","LineWidth",1.5)
end


