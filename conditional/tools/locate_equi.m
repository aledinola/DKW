function jl = locate_equi(xgrid,xi)
% This is the same as locate but assumes the grid is equally spaced
% If xgrid(1)<= xi<=xgrid(N), always returns an index in {1,...,N-1}
% If xi<xgrid(1), then jl=1
% If xi>xgrid(N), then jl=N-1

nx = length(xgrid);

step = xgrid(2)-xgrid(1);
xi_min = xi-xgrid(1);

jl = min(nx-1,max(1,floor(xi_min/step)+1));
 
% if xi<xgrid(1)
%     jl = 0;
%     return
% elseif xi>xgrid(nx)
%     jl = nx;
%     return
% else
%     t = 1 + (xi-xgrid(1))/(xgrid(end)-xgrid(1))*(nx-1);
%     
%     % Matrix element indexing
%     jl=min(max(1,floor(t)),nx-1);
%    
% end


end %END FUNCTION "locate_equi"
