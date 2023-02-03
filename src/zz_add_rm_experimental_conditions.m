% ZZ_ADD_RM_EXPERIMENTAL_CONDITIONS 
%
% inputs:
%   copy_conditions - n x 2 cell array, left column is conditions to copy,
%                       right column is new condition indices (doesn't keep
%                       original copy)
%
%


% Revision history:
%{

2020-03-27: jpm
    Script created

TODO: add everything else


%}


function new_ED = ...
    zz_add_rm_experimental_conditions(Experimental_Data,varargin)

p = inputParser;

% error("Script is broken - need to fix")

% Check for if the entry is a 
isNameOrNameCells = @(x) isstring(x) || ischar(x) || ...
    isstring(x{1}) || ischar(x{1});

addRequired(p,'Experimental_Data',@isstruct);
addParameter(p,'copy_conditions',[],@iscell);
addParameter(p,'ED_fields','',isNameOrNameCells);

p.KeepUnmatched = false;
p.CaseSensitive = false;

parse(p,Experimental_Data,varargin{:})

copy_conditions = p.Results.copy_conditions;
ED_fields = p.Results.ED_fields;

% get list of all fields in struct
all_ED_fields = fieldnames(Experimental_Data);

% IF not given, replace all fields
if isempty(ED_fields)
    ED_fields = all_ED_fields;
elseif isstring(ED_fields) || ischar(ED_fields)
    % Also check if a string was given, and if so, put it into a cell
    ED_fields = {ED_fields};
end

% Get list of fields to copy over as-is
unchanged_fields = setdiff(all_ED_fields,ED_fields);

if ~isempty(copy_conditions)

    % Loop through relevant fields (those that are cells)
    for f = 1:length(ED_fields)
        
        fieldname = ED_fields{f};
        field = Experimental_Data.(fieldname);
        
        if iscell(field)
            
            for row = 1:size(copy_conditions,1)
                
                old_cell = copy_conditions{row,1};
                new_cells = copy_conditions{row,2}(:)';
                
                assert(length(old_cell) == 1,...
                    "The left column of 'copy_conditions' must be integers");
                
                for cell = new_cells
                    
                    new_ED.(fieldname){cell,1} = field{old_cell};
                end
        
            end
            
        else
            % If not a cell, copy over as-is
            new_ED.(fieldname) = field;
            
        end
        
    end
    
    for f = 1:length(unchanged_fields)
        
        new_ED.(unchanged_fields{f}) = ...
            Experimental_Data.(unchanged_fields{f});
    end
        
end

end