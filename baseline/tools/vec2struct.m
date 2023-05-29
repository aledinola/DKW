function [pstruct] = vec2struct(pvec,parnames,pstruct)
%        [par]     = vec2struct(guess,calibNames,par);
%{
 Purpose: adds each element of "pvec" into the structure "pstruct", using 
 the elements in "parnames" as field names.
 Does the same thing as "updatepar" commented below.
 INPUTS:
    "pvec"       Numerical vector
    "parnames"   Character array of strings
    "pstruct"    Structure
 OUTPUTS:
    "pstruct"    Structure (updated)
%}

n = numel(parnames);
for i = 1:n
   pstruct.(parnames{i}) = pvec(i);
end

end %end function "vec2struct"

