function pvec = struct2vec(pstruct,pnames)
    %STRUCT2VEC: Procedure that extract vector from structure with scalar fields. 
    %
    % SYNTAX:
    %   pvec = estim.struct2vec(pstruct,pnames)
    % 
    % OUTPUTS
    %   pvec:      k dimentional vector of scalars
    %
    % INPUTS
    %   pstruct:   Structure with fields taking scalar values (at least one field for each value element in pnames) 
    %
    %   pnames:    k dimensional cell array, that specify names of the fields in pstruct that should be written to pvec
    %
    % See also: estim.vec2struct
	
	%NOTE: To increase speed, comment out these checks on inputs (lines 18-24)
    %NOTE: similar to CreateVectorFromParams2.m written by R.Kirkby
	if ~isstruct(pstruct)
		error('First input must be a structure')
	end	
	
	if ~iscell(pnames)
		error('Second input must be a cell array of strings')
	end
	
	k=numel(pnames);
    j = 1;
    for i=1:k
        r_i = size(pstruct.(char(pnames(i))),1);
        c_i = size(pstruct.(char(pnames(i))),2);
        for c=1:c_i
            for r=1:r_i
                pvec(j,1)=pstruct.(char(pnames(i)))(r, c);
                j = j + 1; 
            end
        end    
    end
end % end of estim.struct2vec
