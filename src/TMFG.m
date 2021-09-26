function [A,tri,separators,cliques,cliqueTree]=TMFG(W)
% TMFG  
%           Computes a Planar Maximally Filtered Graph (PMFG) starting from 
%           a tetrahedron and inserting recursively vertices inside 
%           existing triangles (T2 move) in order to approximate a
%           maxiaml planar graph with the largest total weight - non
%           negative weights
% Function call 
%           [A,tri,clique3]=TMFG(W)
%           [A,tri,clique3,clique4]=TMFG(W)
%           [A,tri,clique3,clique4,cliqueTree]=TMFG(W)
% Input     
%           W:  a NxN matrix of -non-negative- weights 
% Output
%           A: adjacency matrix of the PMFG (with weights)
%           tri: list of triangles (triangular faces)
%           clique3: list of 3-cliques that are not triangular faces (all
%           3-cliques are given by: [tri;clique3])
%           clique4: (optional) list of all 4-cliques
%           cliqueTree: (optional) 4-cliques tree structure (adjacency matrix)
% Example:
%    [A,tri,clique3]=TMFG(corr(randn(80,100)));
%    [A,tri,clique3,clique4]=TMFG(corr(randn(80,100)));
%    [A,tri,clique3,clique4,cliqueTree]=TMFG(corr(randn(80,100)));
%
% Reference: http://arxiv.org/pdf/1505.02445.pdf
% Guido Previde Massara,  T. Di Matteo, and Tomaso Aste. "Network Filtering
% for Big Data: Triangulated Maximally Filtered Graph." 
% arXiv preprint arXiv:1505.02445 (2015).
%
% TA, GPM 10/05/2014
%
N    = size(W,1);
if N< 9, fprintf('W Matrix too small \n'), end
if any(W<0), fprintf('W Matrix has negative elements! \n'), end
A    = sparse(N,N);     % ininzialize adjacency matrix
in_v = zeros(N,1);      % ininzialize list of inserted vertices
tri  = zeros(2*N-4,3);  % ininzialize list of triangles
separators=zeros(N-4,3);   % ininzialize list of 3-cliques (non face-triangles)
%% find 3 vertices with largest strength
s    = sum(W.*(W>mean(W(:))),2);
[~,j]=sort(s,'descend');
in_v(1:4)  = j(1:4);
ou_v = setdiff([1:N],in_v); % list of vertices not inserted yet
%% build the tetrahedron with largest strength
tri(1,:)=in_v([1 2 3]);
tri(2,:)=in_v([2 3 4]);
tri(3,:)=in_v([1 2 4]);
tri(4,:)=in_v([1 3 4]);
A(in_v(1),in_v(2)) = 1; 
A(in_v(1),in_v(3)) = 1;
A(in_v(1),in_v(4)) = 1;
A(in_v(2),in_v(3)) = 1;
A(in_v(2),in_v(4)) = 1;
A(in_v(3),in_v(4)) = 1;
%% build initial gain table
gain = -inf(N,2*N-4);
gain(ou_v,1) = sum(W(ou_v,tri(1,:)),2);
gain(ou_v,2) = sum(W(ou_v,tri(2,:)),2);
gain(ou_v,3) = sum(W(ou_v,tri(3,:)),2);
gain(ou_v,4) = sum(W(ou_v,tri(4,:)),2);
kk = 4;  % number of triangles
for k=5:N
    %% find best vertex to add in a triangle
    if length(ou_v)==1 %special case for the last vertex
        ve = ou_v;
        v  = 1;
        [~,tr] = max(gain(ou_v,:));
    else
        [gij,v]= max(gain(ou_v,:));
        [~,tr] = max( gij );
        ve = ou_v(v(tr));
        v  = v(tr);
    end
    %% update vertex lists
    ou_v = ou_v([1:(v-1),(v+1):end]);
    in_v(k) = ve;
    %% update adjacency matrix
    A(ve,tri(tr,:))=1;
    %% update 3-clique list
    separators(k-4,:) = tri(tr,:); 
    %% update triangle list replacing 1 and adding 2 triangles 
    tri(kk+1,:) = [tri(tr,[1,3]),ve]; % add
    tri(kk+2,:) = [tri(tr,[2,3]),ve]; % add
    tri(tr,:)   = [tri(tr,[1,2]),ve]; % replace
    %% update gain table
    gain(ve,:)=0;
    gain(ou_v,tr)  = sum(W(ou_v,tri(tr,:)),2);
    gain(ou_v,kk+1)= sum(W(ou_v,tri(kk+1,:)),2);
    gain(ou_v,kk+2)= sum(W(ou_v,tri(kk+2,:)),2);
    %% update number of triangles
    kk = kk+2; 
    % if mod(k,1000)==0,fprintf('TMFG - T2 only: %0.2f per-cent done\n',k/N*100);end
end
A = W.*((A+A')==1);
%% computes 4-clique list
if nargout>3     
    cliques = [in_v(1:4)';separators,in_v(5:end)]; 
end
%% computes 4-clique tree (note this may include incomplete cliques!)
if nargout>4
    cliqueTree=sparse(size(cliques,1),size(cliques,1));
    for i=1:size(cliques,1)
        ss = zeros(size(cliques,1),1);
        for k=1:4
            ss=ss+sum((cliques(i,k)==cliques),2);
        end
        cliqueTree(i,ss==3)=1;
    end
end