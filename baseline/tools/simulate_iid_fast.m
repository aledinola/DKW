function [ind_sim, val_sim] = simulate_iid_fast(z_grid,z_prob,dbg)

% Simulate a random draw from a discretized distribution 
%{

INPUTS:
   zGrid: discretized grid for shock z
   zProb: discretized prob vector for z, must be >=0 and sum to 1
  

OUT:
   ind_sim: index for the draw (integer)
   z_sim = zGrid(ind_sim), real value

Last updated by A.Di Nola: 2021-February-18
Adapted from Kindermann's toolbox "get_tomorrow"
%}

%% Input check

n = length(z_grid);

if dbg
   validateattributes(z_prob(:), {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0, '<=', 1, ...
      'size', [n,1]})
   if abs(sum(z_prob) - 1) > 1e-5
      error('Initial prob does not sum to 1');
   end
end

%% Draw random number (this could be done outside the function)
u = rand(1,1);

z_prob_cum = cumsum(z_prob);

ind_sim = locate(z_prob_cum,u)+1;

if (ind_sim<1 || ind_sim>n)
    fprintf('ind_sim = %d \n',ind_sim)
    error('result is out of bounds')
end
        
val_sim = z_grid(ind_sim);

end %END FUNCTION
