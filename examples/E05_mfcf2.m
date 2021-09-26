addpath('../src/')

% 16/05/2019 Guido Previde Massara
% Example usage of the MFCF routine
% See:
% Massara, G. P., & Aste, T. (2019). Learning Clique Forests. 
% arXiv preprint arXiv:1905.02266.

% random symmetric matrix
M = rand(50,50); M = .5 * (M +M');

% initialise the ct_control structure
%  - maximum clique size
ct_control = ct_control_gen(2,2,2);

% this gain function, easy and quick for testing, returns the sum of the
% square of the links between the new node and the separator
gain_function = @gf_sumsquares_gen;

% the MFCF algo
[cliques, separators, peo, parent_cliques] = MFCF2_matrix(M, ct_control, gain_function);

J = LOGO(M, cliques, separators, true);
S = LOGO(M, cliques, separators, false);


L1 = laplacian(inv(M)!=0);
L2 = laplacian(J!=0);

subplot(2,3,1); imagesc(M)
subplot(2,3,2); imagesc(inv(J))
subplot(2,3,3); plot(eigs(L1, 50))
subplot(2,3,4); imagesc(inv(M))
subplot(2,3,5); imagesc(J)
subplot(2,3,6); plot(eigs(L2,50))

