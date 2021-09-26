function [ cliques, separators, peo, parent_cliques ] = MFCF2_matrix( X, ct_control, gain_function )
%MFCF Builds a clique forest from data
% See:
% Massara, G. P., & Aste, T. (2019). Learning Clique Forests. 
% arXiv preprint arXiv:1905.02266.


% Analyse X;

sizeX = size(X);
isSquare = sizeX(1) == sizeX(2);
p = sizeX(2);
max_cliques = p - 1;

% End analyse X;

% Preallocate output;

% cliques is a vector of cliques stored in max_cliques rows and max_clique_size 
% columns
cliques = NaN(max_cliques, ct_control.max_clique_size);

% separators is a vector of separators stored in max_cliques - 1 rows and 
% max clique_size - 1 columns
separators = NaN(max_cliques - 1, ct_control.max_clique_size - 1);

% a perfect elimination order for the graph
peo = [];

% we no longer output the tree as such, but we output the parent clique of 
% every clique 
parent_cliques = NaN(max_cliques, ct_control.max_clique_size);

% End preallocate output;

% init gain table;
GT = init_gain_table(ct_control.max_clique_size);

% init book-keeping
outstanding_nodes = 1:p;
clique_no = 0;
sep_no = 0;

% If matrix is not square make kernel matrix
if ~isSquare 
    error('This version of MFCF2 works only with square (similarity) matrices.')
else
    sums = sum(X.*(X>mean(X(:))),2);
end

% get first ct_control.min_clique nodes as the first clique
[~,j]=sort(sums,'descend');
first_clique = j(1:ct_control.min_clique_size).';
clique_no = 1;
cliques(clique_no, 1:numel(first_clique)) = first_clique;

% the first clique has no parent clique;
% leave the parent clique slot empty;

% remove first clique form outstanding_nodes
outstanding_nodes = outstanding_nodes(~ismember(outstanding_nodes, first_clique));
peo = first_clique.';

% calculate gains for this initial clique
[nodes gains seps] = gain_function(X, first_clique, outstanding_nodes, ct_control);

% how many records ne have to add to the gain table
new_gains = numel(gains);

% calculate the first and last records to update in gain table
from = GT.tot_records + 1;
to   = from + new_gains - 1;

% update count of records in gain table
GT.tot_records = GT.tot_records + new_gains;

idx = (from:to).';
clq_len = numel(first_clique);

% populate cliques in gain table, new_gains replicas of first clique.
GT.cliques(idx, 1:clq_len) = repmat(first_clique,new_gains,1);
if clq_len < ct_control.max_clique_size
    GT.cliques(idx, (clq_len+1):end) = NaN;
end

% populate gains in gain table
GT.gains(idx) = gains;

% populate separators in gain table
GT.separators(idx, :) = seps;

% populate the nodes in gain table
GT.nodes = nodes;

% for all outstanding nodes, consider clique expansions with existing cliques
while(~isempty(outstanding_nodes))

   [the_gain, idx] = nanmax(GT.gains);

   % Clique expansion 
   % ----------------
   
   % case no gain, we must add an isolated clique, let us pick up the first 
   % outstanding node
   if isnan(the_gain)
     the_node = outstanding_nodes(1);
     the_sep  = [];
     parent_clique = [];
     parent_clique_id = NaN;
   else    
     idx = idx(1); % keep only first match
     the_node = GT.nodes(idx);
     the_sep  = GT.separators(idx,:);
     the_sep = the_sep(~isnan(the_sep));
     parent_clique = GT.cliques(idx,:);
     parent_clique_id = id_from_set(cliques, parent_clique);
   end;
   
   new_clique = [the_sep the_node];
   outstanding_nodes = outstanding_nodes(outstanding_nodes ~= the_node);
   peo = [peo; the_node];
   
   %togo = numel(outstanding_nodes)
   
   % track if it is a clique extension or a new clique
   clique_extension = 0;
   
   % in case of no gain add an isolated clique
   if isnan(the_gain)
     clique_no = clique_no + 1;
     clique_to_update = clique_no;
     cliques(clique_to_update, 1:numel(new_clique)) = new_clique;
     % empty parent clique in this case;
   % case new clique with existing intersection
   elseif numel(new_clique) <= sum(~isnan(parent_clique))
       clique_no = clique_no + 1;
       clique_to_update = clique_no;
       sep_no = sep_no + 1;
       cliques(clique_to_update, 1:numel(new_clique)) = new_clique;
       separators(sep_no, 1:numel(the_sep)) = the_sep;
       parent_cliques(clique_to_update, 1:numel(parent_clique)) = parent_clique;
   % case extension of existing clique
   else
       clique_to_update = parent_clique_id;
       cliques(clique_to_update, 1:numel(new_clique)) = new_clique;
       clique_extension = 1;
       % parent clique is the same as the clique without the last estension
       old_clique_idx = id_from_set(GT.cliques, parent_clique);
   end
   
   
   % don't update the gain table when nodes finished
   if isempty(outstanding_nodes) 
       break;
   end;
    
   % finally update gain table
   [nodes gains seps] = gain_function(X, new_clique, outstanding_nodes, ct_control);
   new_gains = numel(gains);
   from = GT.tot_records + 1;
   to   = from + new_gains - 1;
   GT.tot_records = GT.tot_records + new_gains;
   idx = (from:to).';
   new_clique = new_clique(~isnan(new_clique));
   clq_len = numel(new_clique);
   GT.cliques(idx, 1:clq_len) = repmat(new_clique,new_gains,1);
   if clq_len < ct_control.max_clique_size
      GT.cliques(idx, (clq_len+1):end) = NaN;
   end
   GT.gains(idx) = gains;
   GT.separators(idx, :) = seps;
   GT.nodes(idx)= nodes;
   
   % remove from gain table where node is the new one
   idx = find(GT.nodes == the_node);
   GT.gains(idx) = NaN;
   
   % if the clique was expanded remove the records with the old clique from the 
   % gain table
   if clique_extension == 1
     GT.gains(old_clique_idx) = NaN;
   end

   % remove separators who have reached max cooridination number
   if ct_control.coordination_num > 0
       % count appearance of separator just used
       idx = id_from_set(separators, the_sep);
       if numel(idx) >= ct_control.coordination_num
          idx = id_from_set(GT.separators, the_sep);
          GT.gains(idx) = NaN;
       end;
   end
end
   peo = flipud(peo);
end

