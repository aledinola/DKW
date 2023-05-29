function [] = mystruct2table_mom(pstruct1,pstruct2,pnames,calibWeights,targetNames_long,tex,tabDir,filename)

%{
DESCRIPTION:
Generate a table from structures pstruct1 and pstruct2, 
considering only pnames

INPUTS:
- pstruct1:      structure with scalar fields (data_mom)
- pstruct1:      structure with scalar fields (model_mom)
- pnames:        cell array of characters
- calibWeights:  TBA
- targetNames_long: cell array of characters
- tex:           0/1 flag; if 1, create tex file
- tabDir:        Folder where you want to save the latex file
- filename:      Name of the latex file, must have .tex suffix
%}

if ~isstruct(pstruct1)
    error('Input pstruct1 in mystruct2table_mom must be a structure')
end
if ~isstruct(pstruct2)
    error('Input pstruct2 in mystruct2table_mom must be a structure')
end
if ~iscell(pnames)
    error('Input pnames in mystruct2table_mom must be a cell array')
end
if ~iscell(targetNames_long)
    error('Input targetNames_long in mystruct2table_mom must be a cell array')
end
if ~ischar(tabDir)
    error('Input tabDir in mystruct2table_mom must be a character var')
end
if ~ischar(filename)
    error('Input filename in mystruct2table_mom must be a character var')
end

fnames = pnames;

% if length(fieldnames(pstruct1))~=length(fieldnames(pstruct2))
%     error('pstruct1 and pstruct2 must have the same number of fields')
% end

% Check if all fields of pstruct are scalars
for i = 1:numel(fnames)
    aux = size(pstruct1.(fnames{i}));
    if ~isequal(aux,[1,1])
        error("Structure has non-scalar fields, not admissible")
    end
end
for i = 1:numel(fnames)
    aux = size(pstruct2.(fnames{i}));
    if ~isequal(aux,[1,1])
        error("Structure has non-scalar fields, not admissible")
    end
end

% Convert structure into vector
pvec1 = zeros(numel(fnames),1); % data moments
for i = 1:numel(fnames) 
    pvec1(i) = pstruct1.(fnames{i});
end
pvec2 = zeros(numel(fnames),1); % model moments
for i = 1:numel(fnames) 
    pvec2(i) = pstruct2.(fnames{i});
end
calibWeights_vec = zeros(numel(fnames),1); % calibration weights
for i = 1:numel(fnames) 
    calibWeights_vec(i) = calibWeights.(fnames{i});
end

dev = zeros(numel(fnames),1); % abs deviation
for i = 1:numel(fnames) 
    dev(i) = calibWeights_vec(i)*((pvec1(i)-pvec2(i))/pvec1(i))^2;
end


width = max(cellfun('length', fnames));  
width = max([width,length('Moment'),length('Data')]);

%% Display on screen - this we always do
fprintf('--------------------------------------------------------------\n');
fprintf(' %-*s  %-s    %-s    %-s \n',width,'Moment','Data','Model','Sqr Dist');
fprintf('--------------------------------------------------------------\n');
fprintf(' \n');
for i = 1:numel(fnames)
    fprintf('%-*s   %-8.4f %-8.4f %-8.4f \n',width,fnames{i},pvec1(i),pvec2(i),dev(i));
end

%% Write to latex table
if tex==1
    
    FID = fopen(fullfile(tabDir,filename),'w');
    fprintf(FID,' \\begin{tabular}{lcc} \\hline \\hline \n');
    fprintf(FID,' Moment & Data & Model \\\\ \n');
    fprintf(FID,'\\hline \n');
    for i = 1:numel(fnames)
        fprintf(FID,'%s  &  %8.4f & %8.4f \\\\ \n',targetNames_long{i},pvec1(i),pvec2(i));
    end
    fprintf(FID, '\\hline \\hline \n \\end{tabular} \n');
    fclose(FID);
    
end

end %END FUNCTION

