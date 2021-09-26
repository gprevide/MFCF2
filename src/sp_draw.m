function retval = sp_draw (L)
[v, e] = eigs(L);
gplot(L, v(:,[3 2]));
hold on;
gplot(L, v(:,[3 2]),"o");
end
