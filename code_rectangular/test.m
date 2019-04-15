close, clear all
clc
% Driver script for solving the 2D Poisson equation
Globals2D;

% Polynomial order used for approximation 
N = 5;

% Read in Mesh
L= 30;
h= 1;
[nodes, EToV,BCType3, cmdout] = pmesh(h, L);
Nv= size(nodes,1);
K= size(EToV,1);
VX= nodes(:,1)';
VY= nodes(:,2)';

% [Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('circA01.neu');
% BCType(BCType==3)= 6;


tri= EToV;
nodes= [VX',VY'];
K= size(tri,1);
Nn= size(nodes,1);

BCType2= zeros(size(tri));
Nfp= N-1;
edges= [tri(:,[1,2]);tri(:,[2,3]); tri(:,[3,1])];
e1= sort(edges,2);
e0= [e1, repmat((1:K)',3,1),repmat((1:K)',3,1)];
id = (e0(:,1)-1)*Nn + e0(:,2);
e2= sortrows([id,e0],1);
[~,ai,~]= unique(e2(:,1),'rows');
P= 1:size(e2,1);
bi= P(~ismember(P,ai));
e2(bi,5)= e2(bi-1,4);
bind= e2;
bind([bi-1,bi],:)= [];
bind(:,[1,5])= [];
e2(bi-1,:)= [];
n= size(bind,1);
for k=1:n
    v= tri(bind(k,3),:);
    p= sort(find(v==bind(k,1) | v==bind(k,2)))';
    ind= [3 1]*p;
    ind= 1*(ind == 5)+2*(ind == 9)+ 3*(ind==6);
    BCType2(bind(k,3),ind)= 6;
end


a= (BCType2-BCType3);
b= find(a==6);
[I, J]= ind2sub(size(BCType),b);
I
c= EToV(I,:)


figure(1)
hold on
triplot(EToV, VX, VY)
triplot(EToV(bind(:,3),:), VX, VY,'r')
