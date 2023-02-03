% Script to compare equality of shared fields between two structs
% Inputs: two structs
% Outputs: cell of each shared fieldname, and equality

% Revision history:
%{

2020-04-07: jpm
    Script created


%}

function field_equality = compare_structs(struct1,struct2)

struct1_fields = fieldnames(struct1);
field_equality = cell(length(struct1_fields),2);

show_fields = true;

for i=1:length(struct1_fields)
    if isfield(struct2,struct1_fields{i})
        eq_check = isequal(struct1.(struct1_fields{i}),...
                            struct2.(struct1_fields{i}));
        field_equality{i,1} = struct1_fields{i};
        field_equality{i,2} = eq_check;
        
        if show_fields
            field_equality{i,3} = struct1.(struct1_fields{i});
            field_equality{i,4} = struct2.(struct1_fields{i});
        end
        
    end
end

no_match_rows = find(cellfun(@isempty,field_equality(:,2)));

field_equality(no_match_rows,:) = [];

end