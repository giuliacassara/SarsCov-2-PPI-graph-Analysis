m = matfile('sarshuman.mat');
whos('-file','sarshuman.mat');
A = m.sarshumanAdjMat;

G = graph(A,'omitselfloops');
p = plot(G);

L = laplacian(G);
[V,D] = eigs(L,2,'smallestabs');

w = V(:,2);
highlight(p,find(w>=0),'NodeColor','r') % subgraph A
highlight(p,find(w<0),'NodeColor','k') % subgraph B

