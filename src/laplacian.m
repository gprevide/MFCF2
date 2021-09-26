function retval = laplacian (M)
  retval = diag(sum(M)) - M;
end

