function [nodes, tri, BCType, cmdout] = Mesh(h, Lr,Lq,w, g)

points= [-Lq/4-g/2+w/4,Lq/2;-Lq/4-g/2-w/4,Lq/2;-Lq/4-g/2-w/4,-Lq/2;-Lq/4-g/2+w/4,-Lq/2;...  %1-4
    -Lq/4+g/2-w/4,-Lq/2;-Lq/4+w/4+g/2,-Lq/2;Lq/4-g/2-w/4,-Lq/2;Lq/4-g/2+w/4,-Lq/2;...                             %5-8
    Lq/4+g/2-w/4,-Lq/2;Lq/4+g/2+w/4,-Lq/2;Lq/4+g/2+w/4,Lq/2;Lq/4+g/2-w/4,Lq/2;              %9-12
    Lq/4+w/4-g/2,Lq/2;Lq/4-w/4-g/2,Lq/2;-Lq/4+g/2+w/4,Lq/2;-Lq/4+g/2-w/4,Lq/2;...                               %13-16
    -Lr/2,-Lq/2;Lr/2,-Lq/2;Lr/2,Lq/2;-Lr/2,Lq/2;...                                                     %17-20
    0,-g-w/2;-w/4,-g;0,-g+w/2;w/4,-g;...                                                            %21-24
    0,g-w/2;-w/4,g;0,g+w/2;w/4,g;...                                                                %25-28
    g/2,-w/2;g/2-w/4,0;g/2,w/2;g/2+w/4,0;...                                                        %29-32
    -g/2,-w/2;-g/2-w/4,0;-g/2,w/2;-g/2+w/4,0];                                                      %33-36


elements= [2,20;20,17;17,3;3,34;34,21;21,2;3,4;4,33;33,34;33,36;36,35;35,34;...
    35,1;1,2;4,5;5,22;22,34;22,23;23,36;5,6;6,21;21,24;24,23;6,7;7,21;...
    7,8;8,24;8,9;9,29;29,30;30,25;25,26;26,16;16,1;9,10;10,32;32,31;31,28;...
    28,27;27,15;15,16;10,18;18,19;19,11;11,32;11,12;12,31;12,13;13,28;...
    13,14;14,27;14,15]-1;


points= [0:length(points)-1;points'];
elements= [0:length(elements)-1; elements'];
file= sprintf('domain.poly');
fileID= fopen(file,'w');
fprintf(fileID,'#Gebiet\n%d 2 0 0\n',size(points,2));
fprintf(fileID,'%d\t%0.2f\t%0.2f\n', points);
fprintf(fileID,'%d 0\n',size(elements,2));
fprintf(fileID,'%d\t%d\t%d\n',elements);
fprintf(fileID,'0\n');
fclose(fileID);

command= sprintf('triangle -i -C -q -a%f -p domain',h);
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
delete('domain.1.poly')
delete('domain.poly')

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

