function [ JS ] = LOGO( S, cliques, separators, inverse )
%LOGO assembles the direct or inverse (recommended) matrix from cliques and
%separators
N = size(S,1);
JS = sparse(N,N);

for c=cliques.'
    clique = c(~isnan(c));
    if inverse == true
        JS(clique, clique) = JS(clique, clique) + inv(S(clique, clique));
    else
        JS(clique, clique) = JS(clique, clique) + S(clique, clique);
    end
end

for s=separators.'
    separator = s(~isnan(s));
    if inverse == true
        JS(separator, separator) = JS(separator, separator) - inv(S(separator, separator));
    else
        JS(separator, separator) = JS(separator, separator)- S(separator, separator);
    end
end

end



