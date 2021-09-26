addpath('../src/')

% 16/05/2019 Guido Previde Massara
% Example usage of the MFCF routine
% See:
% Massara, G. P., & Aste, T. (2019). Learning Clique Forests. 
% arXiv preprint arXiv:1905.02266.

% random symmetric matrix
M = rand(200,200); M = .5 * (M +M');

% initialise the ct_control structure
%  - maximum clique size
ct_control = ct_control_TMFG();

% this gain function, easy and quick for testing, returns the sum of the
% square of the links between the new node and the separator
gain_function = @gf_sumsquares_gen;

% the MFCF algo
[cliques, separators, peo, tree] = MFCF(M, ct_control, gain_function);

J = LOGO(M, cliques, separators, true);
S = LOGO(M, cliques, separators, false);


subplot(2,2,1); imagesc(M)
subplot(2,2,2); imagesc(inv(J))
subplot(2,2,3); imagesc(inv(M))
subplot(2,2,4); imagesc(J)


