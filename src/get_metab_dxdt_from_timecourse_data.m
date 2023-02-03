% GET_METAB_DXDT_FROM_TIMECOURSE_DATA Return a n x 2 array of the change 
% over time of each metabolite measured in a specified experimental
% condition
% 
%

%{ 
Revision history:
2020-06-12: jpm
    Created script, adapted from adjust_init_metab_concs and
        get_fluxes_from_timecourse_data

2020-09-15
    Also now pass back the metab concs at the specified time
%}
function [dxdt, init_metab_concs, init_metab_inds] = get_metab_dxdt_from_timecourse_data(...
    model,Experimental_Data,Options,varargin)

p = inputParser;
checkBool = @(x) islogical(x) || isnumeric(x);

addRequired(p,'model',@isstruct)
addRequired(p,'Experimental_Data',@isstruct)
addRequired(p,'Options',@isstruct)
addParameter(p,'ref_condition',[],@(x) isnumeric(x) && length(x) == 1)
addParameter(p,'ref_timeframe',[],@(x) isnumeric(x) && length(x) == 2)
addParameter(p, 'set_neg_concs_to_zero', true, checkBool) 

p.KeepUnmatched = false;
p.CaseSensitive = false;

parse(p,model,Experimental_Data,Options,varargin{:})
ref_condition = p.Results.ref_condition;
ref_timeframe = p.Results.ref_timeframe;
set_neg_concs_to_zero = p.Results.set_neg_concs_to_zero;

% Get fields from parser if needed
if isempty(ref_condition)
    ref_condition = Options.ref_from_ED_experimental_condition;
end
if isempty(ref_timeframe)
    ref_timeframe = Options.ref_from_ED_flux_timeframe;
end

% Adapted from adjust_init_metab_concs:
ref_timecourse = Experimental_Data.metab_timecourse{ref_condition};
start_time = ref_timeframe(1);
stop_time = ref_timeframe(2);
measured_mets = ref_timecourse(:,1);

dxdt = zeros(length(measured_mets), 2);
init_metab_concs = zeros(length(measured_mets), 1);
init_metab_inds = zeros(length(measured_mets), 1);

for m = 1:length(measured_mets)
    met_ind = find(strcmpi(model.metabs_and_enzyme_complexes,measured_mets{m}));
    assert(length(met_ind) == 1,"Could not find met_ind")
    
    timepoints = ref_timecourse{m,2}(1,:);
    
    t_start_ind = find(timepoints == start_time);
    t_stop_ind = find(timepoints == stop_time);
    assert(length(t_start_ind) == 1 && length(t_stop_ind) == 1,...
        "Could not find t ind")
    
    % Get indices of timepoints in timecourse
    t_inds = t_start_ind : t_stop_ind;
    metab_concs = ref_timecourse{m,2}(2,t_inds);
    % Unless optioned, make sure all concs are non-negative
    if set_neg_concs_to_zero
        metab_concs = max(zeros(size(metab_concs)), metab_concs);        
    end
    
    timespan = timepoints(t_start_ind : t_stop_ind);
    
    % Set up linear solver to get slope (and intercept)
    t1 = [ones(length(t_inds),1), timespan(:)];
    % y = b1 * t   ==>   b1 = (t1)^-1 * y  = t \ y
    
    b1 = t1 \ metab_concs(:);
    
    % First element of b1 is intercept (times the column of ones), 
    % second is slope
    dxdt(m,2) = b1(2);
    dxdt(m,1) = met_ind;
    
    
    % Also save initial metab concs
    init_metab_concs(m) = metab_concs(1);
    init_metab_inds(m) = met_ind;
    
end

end