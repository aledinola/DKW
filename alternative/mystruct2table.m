function [] = mystruct2table(pstruct,pnames,description,dispNames,header,tex,tabDir,filename)
%{
DESCRIPTION:
Generate a table from structure pstruct, considering only pnames
(if pnames is present)

INPUTS:
- pstruct: structure with scalar fields
- pnames: cell array of characters, with names of fields to display
- description: cell array of characters, with more information
- dispNames: cell array of characters
- header: title for columns
- tex: 0/1 flag; if 1, create tex file
- tabDir: folder where to save the tex table (if tex=1)
- filename: name of the tex file
%}

if ~iscell(pnames)
    error('Input pnames in mystruct2table must be a cell array of characters')
end
if ~iscell(description)
    error('Input description in mystruct2table must be a cell array of characters')
end
if ~iscell(dispNames)
    error('Input dispNames in mystruct2table must be a cell array of characters')
end

if ~ischar(tabDir)
    error('Input tabDir in mystruct2table must be a character var')
end
if ~ischar(filename)
    error('Input filename in mystruct2table must be a character var')
end

% If there is only one input, use all the field names of pstruct
if isempty(pnames)
    fnames = fieldnames(pstruct);
else
    fnames = pnames;
end

% Check that character arrays pnames and description (if present) have the 
% same number of elements
if ~isempty(description)
    if ~isequal(numel(fnames),numel(description))
        error("character arrays pnames and description MUST have the same number of elements")
    end
end

% Check if all fields of pstruct are scalars
for i = 1:numel(fnames)
    aux = size(pstruct.(fnames{i}));
    if ~isequal(aux,[1,1])
        error("Structure has non-scalar fields, not admissible")
    end
end

% Convert structure into vector
pvec = zeros(numel(fnames),1);
for i = 1:numel(fnames) 
    pvec(i) = pstruct.(fnames{i});
end

width = max(cellfun('length', fnames));
if ~isempty(header)
    width = max(width, length(header{1}));
    fprintf('------------------------------\n')
    fprintf('%-*s   %s \n', width, header{1}, header{2});
    fprintf('------------------------------\n')
end
% Display table with parameter i name in string fnames{i} and numerical value in
% pvec(i)
for i = 1:numel(fnames)
    fprintf('%-*s   %-8.16f\n', width, fnames{i},pvec(i));
end

if tex==1
    
    
    
    FID = fopen(fullfile(tabDir,filename),'w');
    
    if isempty(description)
        % Without extra column with description of parameters
        fprintf(FID,' \\begin{tabular}{lc} \\hline \\hline \n');
        fprintf(FID,' Parameter & Description & Value \\\\ \n');
        fprintf(FID,' \\hline \n');
        for i = 1:numel(fnames)
            fprintf(FID,'%s  &  %8.3f \\\\ \n',dispNames{i},pvec(i));
        end
        fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
        
    else
        % With extra column with description of parameters
        fprintf(FID,' \\begin{tabular}{llc} \\hline \\hline \n');
        fprintf(FID,' Parameter & Description & Value \\\\ \n');
        fprintf(FID,' \\hline \n');
        for i = 1:numel(fnames)
            fprintf(FID,'%s  & %s & %8.3f \\\\ \n',dispNames{i},description{i},pvec(i));
        end
        fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
        
    end
    
    fclose(FID);
    
    
    
end

end %END FUNCTION

