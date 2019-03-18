function [nodes, tri, BCType, cmdout] = pmesh(h,L)

points= [-L/2,-L/2;L/2,-L/2;L/2,L/2;-L/2,L/2];
points= [0:length(points)-1;points'];
file= sprintf('domain.node');
fileID= fopen(file,'w');
fprintf(fileID,'#Gebiet\n%d 2 0 0\n',size(points,2));
fprintf(fileID,'%d\t%0.2f\t%0.2f\n', points);
fclose(fileID);

command= sprintf('triangle -i -C -q -a%f domain',h);
[status,cmdout] = system(command);

%% Lese Mesh aus

file= sprintf('domain.1.node');
fileID= fopen(file,'r');
first= fscanf(fileID,'%d',[1 4]);
second= textscan(fileID, '%f');
skip= sum(first(2:4))+1;
nodes= [second{1}(2:skip:end), second{1}(3:skip:end)];
marker= second{1}(4:skip:end);
fclose(fileID);
delete(file)
file= sprintf('domain.1.ele');
fileID= fopen(file,'r');
first= fscanf(fileID,'%d',[1 3]);
second= textscan(fileID, '%f');
skip= sum(first(2:3))+1;
tri= [second{1}(2:skip:end), second{1}(3:skip:end), second{1}(4:skip:end)]+1;
fclose(fileID);
delete(file)
delete('domain.node')



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

