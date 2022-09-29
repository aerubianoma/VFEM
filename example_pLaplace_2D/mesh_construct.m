function mesh = mesh_construct(n,plot_mesh)

nodes2coord = [0 0; 1 0; 1 1; 0 1];   % square mesh coordinates
elems2nodes = [1 2 3 ; 1 3 4 ];       % square mesh elements
%edgesNeumann2nodes = [];              % no Neummann boundary conditions => full Dirichlet boundary conditions

edgesDirichlet2nodes=[4 1; 2 3]; edgesNeumann2nodes = [1 2; 3 4];   %two oposite edges
%edgesDirichlet2nodes=[2 3; 3 4]; edgesNeumann2nodes = [4 1; 1 2] ;   %two oposite edges
edgesDirichlet2nodes=[1 2; 2 3; 3 4; 4 1];  edgesNeumann2nodes = [];

%mesh L-shape
%nodes2coord=[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2];
%elems2nodes=[1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7];  
%edgesDirichlet2nodes=[1 2; 2 3; 3 4; 4 5; 5 8; 8 7; 7 6; 6 1]; edgesNeumann2nodes = [];

%uniform mesh refinement
for r=1:n
    [nodes2coord,elems2nodes,edgesDirichlet2nodes,edgesNeumann2nodes] = refinement_uniform_2D(nodes2coord,elems2nodes,edgesDirichlet2nodes,edgesNeumann2nodes);
end
[edges2nodes, edges2elems, ~, nodes2edges] = getEdges(elems2nodes);
nodes2elems = entryInWhichRows(elems2nodes);

% ned = size(edges2nodes,1);
nn = size(nodes2coord,1);  % number of nodes
ne = size(elems2nodes,1);  % number of elements
dim = size(nodes2coord,2); % problem dimension


%all gradients
%particular part for a given element in a given dimension
NLB = 3; %number of local basic functions, it must be known!
coord = zeros(dim,NLB,ne);
for d=1:dim
    for i=1:NLB
        coord(d,i,:) = nodes2coord(elems2nodes(:,i),d);
    end
end   
IP = [1/3;1/3]; [dphi,jac] = phider(coord,IP,'P1'); 
dphi = squeeze(dphi); %all gradients 
areas = abs(squeeze(jac))/factorial(dim); %all areas

%testing gradient
Grad1_elems = squeeze(dphi(1,:,:))';
Grad2_elems = squeeze(dphi(2,:,:))';

if plot_mesh
    figure; 
    show_mesh(elems2nodes,nodes2coord,unique(edgesDirichlet2nodes)); 
    xlabel('x'); ylabel('y'); axis image;
end

edgesBoundary2edges = (edges2elems(:,2)==0);           % logical - all boundary edges
edgesBoundary2nodes = edges2nodes(edgesBoundary2edges,:);

%===================================================================

mesh.dim=dim;
mesh.nodes2coord = nodes2coord;
mesh.nodes2elems = nodes2elems;
mesh.elems2nodes = elems2nodes;
% mesh.nodes2nodes = nodes2nodes;
mesh.nodesBoundary = unique(edgesBoundary2nodes);
mesh.nodesInternal = setdiff((1:nn)',mesh.nodesBoundary);
mesh.nodesFree = setdiff((1:nn)',unique(edgesDirichlet2nodes));

%mesh.nodesFree=
mesh.Grad1_elems = Grad1_elems;
mesh.Grad2_elems = Grad2_elems;
mesh.areas = areas;

%===================================================================

ni = numel(mesh.nodesInternal);

nodesInternal_inv = zeros(nn,1);
for i=1:ni
    j = mesh.nodesInternal(i);
    nodesInternal_inv(j) = i;
end

Hstr=sparse(edges2nodes(:,1),edges2nodes(:,2),ones(size(edges2nodes(:,1),1),1),nn,nn) ...
     +sparse(edges2nodes(:,2),edges2nodes(:,1),ones(size(edges2nodes(:,1),1),1),nn,nn) ...
     +sparse(1:nn,1:nn,ones(nn,1),nn,nn);
%Hstr=Hstr(mesh.nodesInternal,mesh.nodesInternal); 
mesh.Hstr = Hstr;
