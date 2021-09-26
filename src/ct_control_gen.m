% Author: Guido Previde Massara
% Date: 26 Sept 2021
% Returns a ct_control structure that reproduces the TMFG

function retval = ct_control_gen (max_cl, min_cl, coord_num)
% max clique size = 4 
retval.max_clique_size = max_cl;
% min clique size = 4
retval.min_clique_size = min_cl;
% edges with score above threshold will appear in the final result
% no threshold for TMFG-like
retval.threshold = 0.0;
%  cooordination number (use each separator only once)
retval.coordination_num = coord_num;
end