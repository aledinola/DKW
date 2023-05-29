function jl = locate(xx,x)
%function jl = locate(xx,x)
%
% x is between xx(jl) and xx(jl+1)
%
% jl is always in {1,2,..,n-1}
%
% xx is assumed to be monotone increasing

n = length(xx);
if x<xx(1)
    jl = 0;
elseif x>xx(n)
    jl = n;
else
    jl = 1;
    ju = n;
    while (ju-jl>1)
        jm = floor((ju+jl)/2);
        if x>=xx(jm)
            jl = jm;
        else
            ju=jm;
        end
    end
end

jl = min(n-1,max(1,jl));



