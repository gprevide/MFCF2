function [M] = clique_matrix(cliques, T)
to = max(horzcat(cliques{:}));
X = randn(T, to);
    for c = 1:numel(cliques)
        nodes = cliques{c};
        %X(1:T, nodes) = X(1:T, nodes) + (1+randn(T, 1)) * randn(1, numel(nodes));
        X(1:T, nodes) =  X(1:T, nodes) + (1+randn(T, 1)) * ones(1, numel(nodes));
    end
M = corr(X);
end
