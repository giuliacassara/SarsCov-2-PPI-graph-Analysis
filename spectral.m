% GET ADJACENCY MATRIX of SARS-HUMAN-drug graph
clc
clear all

addpath(genpath(pwd))
pwd


%
L = laplacian(g);

%Histogram of eigenvalue spacings.
spectral_partition(g);
plot_eigenvalues(L, 100);
import_dvh_interactome();
h2 = WattsStrogatz(6387,18,0.006);
plot_interactome(h2)

[erdosAdj,n,m] = create_ER_Graph(6387,0.0056)
er = graph(erdosAdj); 
plot_interactome(er)

g = importNet('./dataset/sars_human_human_drug.txt', false);
g = rmedge(g, 1:numnodes(g), 1:numnodes(g)); %remove self-loops from graph
adjMatrix = adjacency(g);

n=1;
fid = fopen('./dataset/node_names.txt','rt'); % 'rt' means "read text"
ca = cell(1, 6387);
while n < 6388
      line = fgetl(fid);
      ca{n} = line;
      n = n + 1;
      %ca{n} = n;
end
fclose(fid);

g.Nodes.Name = ca';
g.Nodes.Degree = centrality(g,'degree');
g.Nodes.Betweenness = centrality(g,'betweenness');
g.Nodes.Closeness = centrality(g,'closeness');
g.Nodes.Pagerank = centrality(g,'pagerank');
plot_interactome(g)

[erdosAdj,n,m] = create_ER_Graph(6387,0.0056)
er = graph(erdosAdj);  
plot_interactome(er)
plot_centrality_graph(er)
%now build a watts strogatz graph
plot_ws_distributions()
plot_erdos_renyi_distributions()

ba = scale_free(6387, 5, 4)
%Plot the $\beta = 0.15$ Watts-Strogatz model graph, 
    %making the size and color of each node proportional to its degree. 

plot_barabasi_distribution(ba)

function plot_interactome(g)

    figure;
    deg = degree(g);
    nSizes = 2*sqrt(deg-min(deg)+0.2);
    nColors = deg;
    p2 = plot(g,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.001)  
    layout(p2,'force','UseGravity',true)
    axis equal
    %title('Erdos Renyi Random graph $N = 6387$ nodes','Interpreter','latex')
    colorbar

%legend('Sars-CoV-2', 'Human', 'Drug')
    L = laplacian(g);
    [V,D] = eigs(L,500); 
    eigenvalues = diag(D);
    [eigenvalues, ind] = sort(eigenvalues, 'ascend');
    eigenvalues = eigenvalues(1:500)
    %plot_eigenvalues(L, 500)
    %plot_distribution(g)
    %plot_histfit(g, 'poisson', eigenvalues)
    %plot_histfit(g, 'weibull', eigenvalues)
    %plot_histfit(ba, 'beta', eigenvalues)
    %plot_histfit(ba, 'gamma', eigenvalues)

    spectral_partition(g)
end


function plot_eigenvalues(L, limiteigs)
    
    [V,D] = eigs(L,limiteigs); 
    eigenvalues = diag(D);
    [eigenvalues, ind] = sort(eigenvalues, 'ascend');
    eigenvalues = eigenvalues(1:limiteigs)
    eigenvalues
    figure;
    hold on;
    h = histogram(eigenvalues)
    legend(sprintf('Eigvals'))
    ylabel('Counts');
    xlabel('Eigenvalues')
    h.Normalization = 'countdensity';
    h.NumBins = 50;
    
    
end

function spectral_partition(g)
    figure
    p = plot(g);
    set(p,'EdgeAlpha',0.01);
    layout(p,'force','UseGravity',true)
    axis equal
    L = laplacian(g);
    [V,D] = eigs(L,2,'smallestabs');
    
    %find the fiedler vector
    w = V(:,2);
    %partition the graph  
    highlight(p,find(w>=0),'NodeColor','r') % subgraph A
    highlight(p,find(w<0),'NodeColor','k') % subgraph B
    
    % sign cut does not always produce a balanced cut.
    % It is always possible to bisect a graph by calculating the median of 
    %w and using it as a threshold value. This partition is called the median cut, 
    %and it guarantees an equal number of nodes in each subgraph.
    w_med = w - median(w);
    figure
    p1 = plot(g);
    set(p1,'EdgeAlpha',0.01);
    layout(p1,'force','UseGravity',true)
    axis equal
    highlight(p1,find(w>=0),'NodeColor','r') % subgraph A
    highlight(p1,find(w<0),'NodeColor','k') % subgraph B

end


function plot_distribution(g)
    %degree
    figure
    hold on
    [~,edges] = histcounts(log10(degree(g)));
    histogram(degree(g),10.^edges)
    set(gca, 'xscale','log')
    %histogram(degree(g),'BinMethod','integers','FaceAlpha',0.9);
    hold off
    %title('Node degree distributions for Drug-Virus-Host interactome')
    xlabel('Degree Centrality Score (log10 scale)')
    ylabel('Number of nodes')
    
    %betweenness distribution 
    
    figure
    hold on
    [~,edges] = histcounts(log10(centrality(g,'betweenness')));
    histogram(centrality(g,'betweenness'),10.^edges)
    set(gca, 'xscale','log')
    %histogram(degree(g),'BinMethod','integers','FaceAlpha',0.9);
    hold off
    %title('Node degree distributions for Drug-Virus-Host interactome')
    xlabel('Betweenness Centrality score (log10 scale)')
    ylabel('Number of nodes')
    
    %Closeness distribution 
    
    figure
    hold on
    [~,edges] = histcounts(log10(centrality(g,'closeness')));
    histogram(centrality(g,'closeness'),10.^edges)
    set(gca, 'xscale','log')
    %histogram(degree(g),'BinMethod','integers','FaceAlpha',0.9);
    hold off
    %title('Node degree distributions for Drug-Virus-Host interactome')
    xlabel('Closeness Centrality Score (log10 scale)' )
    ylabel('Number of nodes')
    
    %Eigenvector distribution 
    
    figure
    hold on
    [~,edges] = histcounts(log10(centrality(g,'eigenvector')));
    histogram(centrality(g,'eigenvector'),10.^edges)
    set(gca, 'xscale','log')
    %histogram(degree(g),'BinMethod','integers','FaceAlpha',0.9);
    hold off
    %title('Node degree distributions for Drug-Virus-Host interactome')
    xlabel('Eigenvector Centrality Score (log10 scale)' )
    ylabel('Number of nodes')
    
    figure
    hold on
    [~,edges] = histcounts(log10(centrality(g,'pagerank')));
    histogram(centrality(g,'pagerank'),10.^edges)
    set(gca, 'xscale','log')
    %histogram(degree(g),'BinMethod','integers','FaceAlpha',0.9);
    hold off
    %title('Node degree distributions for Drug-Virus-Host interactome')
    xlabel('Pagerank Centrality Score (log10 scale)' )
    ylabel('Number of nodes')
    
    %Hubs distribution 
   
end

function plot_histfit(g, method, eigenvalues)
    
    figure;
    hold on;
    %h = histogram(eigenvalues, 'Normalization', 'pdf')
    histfit(eigenvalues, 50, method)
    legend('Eigenvalue distribution fitted with', sprintf(method))
    ylabel('Counts');
    xlabel('Eigenvalues')
    %h.Normalization = 'countdensity';
    
end

function plot_barabasi_distribution(ba)

    figure;
    colormap hsv
    deg = degree(ba);
    nSizes = 2*sqrt(deg-min(deg)+0.2);
    nColors = deg;
    p2 = plot(ba,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.01)      
    layout(p2,'force','UseGravity',true)
    axis equal
    title('Barabasi Model scale-free graph with $N = 6387$ nodes, $m0 = 5$, and $\m = 4$','Interpreter','latex')
    colorbar

    plot_distribution(ba);
    L = laplacian(ba);
    [V,D] = eigs(L,500); 
    eigenvalues = diag(D);
    [eigenvalues, ind] = sort(eigenvalues, 'ascend');
    eigenvalues = eigenvalues(1:500)
    plot_eigenvalues(L, 500)
    plot_histfit(ba, 'poisson', eigenvalues)
    plot_histfit(ba, 'weibull', eigenvalues)
    plot_histfit(ba, 'beta', eigenvalues)
    plot_histfit(ba, 'gamma', eigenvalues)

    spectral_partition(ba)
end

function plot_ws_distributions()
    h2 = WattsStrogatz(6387,18,0.15);
    %plot_distribution(h2);
    L = laplacian(h2);
    [V,D] = eigs(L,500); 
    eigenvalues = diag(D);
    [eigenvalues, ind] = sort(eigenvalues, 'ascend');
    eigenvalues = eigenvalues(1:500)
    spectral_partition(h2)

    %plot_eigenvalues(L, 500)
    %plot_histfit(h2, 'poisson', eigenvalues)
    %plot_histfit(h2, 'weibull', eigenvalues)
    %plot_histfit(h2, 'beta', eigenvalues)
    %plot_histfit(h2, 'gamma', eigenvalues)

end


function plot_ws()
   
    h2 = WattsStrogatz(6387,18,0.15);
    plot(h2,'NodeColor','k','EdgeAlpha',0.01)
    title('Watts-Strogatz Graph with $N = 6387$ nodes, $K = 18$, and $\beta = 0.15$', ...
        'Interpreter','latex')

    h3 = WattsStrogatz(6387,18,0.50);
    plot(h3,'NodeColor','k','EdgeAlpha',0.01)
    title('Watts-Strogatz Graph with $N = 6387$ nodes, $K = 18$, and $\beta = 0.50$', ...
        'Interpreter','latex')
    
    h4 = WattsStrogatz(6387,18,1);
    plot(h4,'NodeColor','k','EdgeAlpha',0.01)
    title('Watts-Strogatz Graph with $N = 6387$ nodes, $K = 18$, and $\beta = 1$', ...
        'Interpreter','latex')
    
    histogram(degree(h2),'BinMethod','integers','FaceAlpha',0.9);
    hold on

    histogram(degree(h3),'BinMethod','integers','FaceAlpha',0.9);
    histogram(degree(h4),'BinMethod','integers','FaceAlpha',0.8);
    hold off
    title('Node degree distributions for Watts-Strogatz Model Graphs')
    xlabel('Degree of node')
    ylabel('Number of nodes')
    legend('\beta = 1.0','\beta = 0.50','\beta = 0.15','Location','NorthWest')
    
    
    %hub formation
    n = 500;
    d = [
         mean(mean(distances(h2))), nnz(degree(h2)>=n);
         mean(mean(distances(h3))), nnz(degree(h3)>=n);
         mean(mean(distances(h4))), nnz(degree(h4)>=n)];
    T = table([0.15 0.50 1]', d(:,1), d(:,2),...
    'VariableNames',{'Beta','AvgPathLength','NumberOfHubs'})

    %Plot the $\beta = 0.15$ Watts-Strogatz model graph, 
    %making the size and color of each node proportional to its degree. 
    colormap hsv
    deg = degree(h4);
    nSizes = 2*sqrt(deg-min(deg)+0.2);
    nColors = deg;
    p2 = plot(h4,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.01)      
    layout(p2,'force','UseGravity',true)
    axis equal
    title('Watts-Strogatz Graph with $N = 6387$ nodes, $K = 18$, and $\beta = 1$','Interpreter','latex')
    colorbar
    
end

function h = WattsStrogatz(N,K,beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end

 
function [G,n,m] = create_ER_Graph(n,p)
%% 
%    Description:
%        this function create Erdos-Renyi random Graph*

%    Output Arguments:
%        G : generated random graph
%        n : graph size, number of vertexes, |V|
%        m : graph size, number of edges, |E|
%    Input Arguments:
%        n : graph size, number of vertexes, |V|
%        p : the probability p of the second definition of Erdos-Renyi model.
seed = n;
rng(seed);
G = spones(triu(sprand(n,n,p),1));
m = nnz(G);
G = G + G';

end

function BA = scale_free(n, m0, m)
% Use the Barabasi-Albert model to generate a scale free graph of size n (as
% described in Albert-Laszlo Barabasi & Reka Albert: "Emergence of scaling
% in random networks")
%
% INPUT
% n: [1]: number of nodes
% m0: [1]: number of initially placed nodes
% m: [1]: number of nodes a new added node is connected to, 1 <= m < m0
%
% OUPUT
% A: [n n] sparse symmetric adjacency matrix representing the generated graph

% Start with a graph of size m0 and add edges to this graph. Each of these m0
% nodes is connected to at least m nodes.
B = zeros(m0, m0);
for i = 1:m0
    neighbors = randsample(m0-1, m);
    neighbors = neighbors + (neighbors>=i);
    B(i,neighbors) = 1;
    B(neighbors,i) = 1;
end

% Create a vector of edges added so far, i.e. nodes edge(2*i) and edge(2*i-1),
% 1 <= i <= nEdges, are connected by an edge.
[rows, columns] = find(triu(B));
nEdges = size(rows, 1);
edges = reshape([rows';columns'], 2*nEdges, 1);
edges = [edges; zeros(2*(n-m0)*m,1)];

% Add nodes m0+1:n, one at a time. Each node is connected to m existing nodes, 
% where each of the existing nodes is chosen with a probability that is
% proportional to the number of nodes it is already connected to. 
used = zeros(n, 1); % is a node already used in a timestep?
for i = m0+1:n
    neighbors = zeros(1, m);
    for j=1:m
       k = edges(randi(2*nEdges));
       while used(k)
           k = edges(randi(2*nEdges));
       end
       used(k) = 1;
       neighbors(j) = k;
    end
    used(neighbors) = 0;
    edges(2*nEdges+1:2*nEdges+2*m) = reshape([repmat(i, 1, m); neighbors], ...
     1, 2*m);
    nEdges = nEdges+m;
end

% finally construct a symmetric adjacency matrix using the vector of edges
edges = edges(1:2*nEdges);
first = edges(1:2:end);
second = edges(2:2:end);
A = sparse([first;second], [second;first], ones(2*nEdges, 1), n, n);

BA = graph(A,'upper')
end % scale_free(...)


