% FIND_BEST_KSETS Find the ksets with the lowest/best fitness scores, and
% return their indices and fitnesses.
% 

% function [val_lowest_max, index_lowest_max, val_lowest_avg, index_lowest_avg]...
%     = find_best_Ksets(Perturbation_Results, top_n_best)
%
% function [val_lowest_max, index_lowest_max, val_lowest_avg, index_lowest_avg]...
%     = find_best_Ksets(fitnesses, top_n_best)

% Revision History:
%{

2019-07-12: jpm
    Created revision history
    Removed normalization for fitnesses in each condition - was really
        messing things up
2019-10-07
    Added parser, removed profanity
2020-08-07
    Removed ability to change fitness fxn within this script

TODO:
- should really have input parser to deal with either PR struct or
    fitnesses being passed in

%}

function [val_lowest_max, index_lowest_max, val_lowest_avg, index_lowest_avg]...
    = find_best_Ksets(PR_or_fits, top_n_best, varargin)

%% Set up parser
p = inputParser;

% make a check function for either string or char array
checkName = @(x) ischar(x) || isstring(x);
checkArrayOrAll = @(x) isnumeric(x) || strcmpi(x,'all');
checkBool = @(x) islogical(x) || x==1 || x==0;

% Add each argument into the parser
addRequired(p,'PR_or_fits',@(x) isstruct(x) || isnumeric(x));
addRequired(p,'top_n_best',@isnumeric)
addParameter(p,'Options',[],@isstruct)
addParameter(p,'fitness_fxn',[],checkName)

p.KeepUnmatched = false;
p.CaseSensitive = false;

parse(p,PR_or_fits,top_n_best,varargin{:})

% Get vars from p.Results struct
Options = p.Results.Options;
fitness_fxn = p.Results.fitness_fxn;

% If fitness_fxn is given, require PR, get new fitness values
if ~isempty(fitness_fxn)
    
    % Throw error - recalculating fitnesses and getting new PR should be
    % done before calling this
    error("Recalculating fitnesses and getting new PR should be done before calling `find_best_ksets`\n")
   
    % Need PR_or_fits to be PR struct
    if ~isstruct(PR_or_fits)
        error("'PR_or_fits' field must be PR to recalc fitnesses with new fxn");
    end
    % Also going to need Options
    if isempty(Options)
        error("Must have 'Options' as varargin to recalc fitness with new fxn");
    end
    
    PR_or_fits = POST_recalc_kset_fitness(model,Experimental_Data,...
            PR_or_fits,Options,'fitness_fxn',fitness_fxn);
    
end

% check if Perturbation_Results or fitnesses was passed in
if isstruct(PR_or_fits)
    if isfield(PR_or_fits,'Ksets_kept')
        n_ksets_remaining = length(PR_or_fits.Ksets_kept{end});
    else
        n_ksets_remaining = length(PR_or_fits.avg_fitness_values);
    end
    fit = PR_or_fits.fitness_values;
else
    n_ksets_remaining = size(PR_or_fits,1);
    fit = PR_or_fits;
end

% Get indices of best
largest_per_Kset = (max(fit,[],2));
avg_per_Kset = mean(fit,2, 'omitnan');

[val_lowest_max, index_lowest_max] = sort(largest_per_Kset);
[val_lowest_avg, index_lowest_avg] = sort(avg_per_Kset);

if (n_ksets_remaining >= top_n_best)
    
    val_lowest_max = val_lowest_max(1:top_n_best);
    index_lowest_max = index_lowest_max(1:top_n_best);
    val_lowest_avg = val_lowest_avg(1:top_n_best);
    index_lowest_avg = index_lowest_avg(1:top_n_best);
    
end


end
