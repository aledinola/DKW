function [ind_sim, val_sim] = simulate_iid(z_grid,z_prob,u,dbg)

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
%u = rand(1,1);

for i = 1:n-1
    %Cumulative sum p1+p2+..+pi
    prob_sum = sum(z_prob(1:i));
    if (u <= prob_sum )
        ind_sim = i;
        val_sim = z_grid(ind_sim);
        return
    end %END IF
end %end for

% Else, choose the last value
ind_sim = n;
val_sim = z_grid(ind_sim);

end %END FUNCTION
