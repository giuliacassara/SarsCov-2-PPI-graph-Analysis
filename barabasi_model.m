
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code adopted from                                     %%%
%%% Modeling and Simulating Social Systems with MATLAB    %%%
%%% http://www.soms.ethz.ch/teaching/MatlabFall2012       %%%
%%% Authors: Stefan Brugger and Cristoph Schwirzer, 2011  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adjscalefree = scalefree(6387, 5, 5)



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
    
    figure;
    hold on;
    %h = histogram(eigenvalues, 'Normalization', 'pdf')
    histfit(eigenvalues, 50, 'weibull')
    %legend(sprintf('Eigvals'))
    ylabel('Counts');
    xlabel('Eigenvalues')
    %h.Normalization = 'countdensity';
    
end
