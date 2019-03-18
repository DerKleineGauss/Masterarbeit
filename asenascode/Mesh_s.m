function [nodes, tri, BCType] = Mesh_s(Lr,Lq, ncellx, ncelly)

x= linspace(-Lr/2,Lr/2,ncellx);
y= linspace(-Lq/2,Lq/2,ncelly);
[X, Y]= meshgrid(x,y);
tri= delaunay(X,Y);
nodes= [X(:), Y(:)];

%% BCType
K= size(tri,1);
Nn= size(nodes,1);
BCType= zeros(size(tri));
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
n= size(bind,1);
for k=1:n
    v= tri(bind(k,3),:);
    p= sort(find(v==bind(k,1) | v==bind(k,2)))';
    ind= [3 1]*p;
    ind= 1*(ind == 5)+2*(ind == 9)+ 3*(ind==6);
    BCType(bind(k,3),ind)= 6;
end



end

