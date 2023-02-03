% FUNCTION CHECK_FITNESS_FUNCTION Makes sure that the specified fitness
% function has all requisite info, and prompts user if not
%
%

% Revision history:
%{
2019-09-18: jpm
    Script created
2019-10-07
    Parser to allow either Options or just fitness_fxn name

%}

function Experimental_Data = ...
    check_fitness_function(model,Experimental_Data,varargin)

p = inputParser;
checkName = @(x) ischar(x) || isstring(x);

addRequired(p,'model',@isstruct)
addRequired(p,'Experimental_Data',@isstruct)
addParameter(p,'Options',[],@isstruct)
addParameter(p,'fitness_fxn','',checkName)

p.KeepUnmatched = false;
p.CaseSensitive = false;

parse(p,model,Experimental_Data,varargin{:})

Options = p.Results.Options;
fitness_fxn = p.Results.fitness_fxn;

% Only use fitness_fxn if Options is empty
if isempty(Options)
    if isempty(fitness_fxn)
        error("Need either 'Options' or 'fitness_fxn' as varargin");
    end 
        
else
    fitness_fxn = Options.fitness_fxn;
end

% Get indices of all reactions used in fitness fxn across all conditions
n_conditions = length(Experimental_Data.flux_rxns);
all_exp_rxns = [];
for cond = 1:n_conditions
    all_exp_rxns = [all_exp_rxns; Experimental_Data.flux_rxns{cond}(:)];  
end
all_exp_rxns = unique(all_exp_rxns);

if strcmpi(fitness_fxn,'wgt_exp_error')
    % 'rxn_weights' should have an element for each rxn
    if isfield(Experimental_Data,'rxn_weights') && ...
            (length(Experimental_Data.rxn_weights) == length(model.rxns) ||...
                length(Experimental_Data.rxn_weights) == size(model.S,2))
        
        rxn_weights = Experimental_Data.rxn_weights(:);
    else
        rxn_weights = zeros(length(model.rxns),1);
        
        input_rxn_inds = input('Input an array of reactions (by index) with non-zero weights:'); 
        input_rxn_weights = input('Input the weights of those rxns (SAME LENGTH): ');
        
        if length(input_rxn_inds) ~= length(input_rxn_weights)
            error("Inputted weights are not the same length as the number of rxns");
        end
        
        rxn_weights(input_rxn_inds) = input_rxn_weights;
                
    end
    
    Experimental_Data.rxn_weights = rxn_weights;
        
    
elseif strcmpi(fitness_fxn,'wgt_rmse')
    error("TODO: not implemented")    
    
end


end