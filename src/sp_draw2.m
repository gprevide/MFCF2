function retval = sp_draw2 (L)
[v, e] = eigs(L);
gplot(L, v(:,[4 3]));
hold on;
gplot(L, v(:,[4 3]),"o");
end
