% Function to plot a confidence interval, given timecourses for several
% ksets for a single run/condition


function [h_CI, h_best] = plot_timecourse_CI_intervals(tc_x_vals, all_kset_tc, ...
    kset_inds, base_color,...
    varargin)

p = inputParser;

% Set up input parser
isName = @(x) isstring(x) || ischar(x) || iscell(x);
isBool = @(x) isnumeric(x) || islogical(x);

CI_int_default = [0.95, 0.5];

addRequired(p, 'tc_x_vals', @isnumeric)
addRequired(p, 'all_ksets_tc', @isnumeric)
addRequired(p, 'kset_inds', @isnumeric)
addRequired(p, 'base_color', @isnumeric)
addParameter(p, 'CI_int', [], @isnumeric)
addParameter(p, 'plot_best_kset', true, isBool)
addParameter(p, 'condition_name', [], isName)

p.KeepUnmatched = false;
p.CaseSensitive = false;

% Allow inputing all varargin in a single cell
if length(varargin) == 1
    varargin = varargin{1};
end

parse(p, tc_x_vals, all_kset_tc, kset_inds, base_color, varargin{:});

CI_int = p.Results.CI_int;
condition_name = p.Results.condition_name;
plot_best_kset = p.Results.plot_best_kset;

if isempty(CI_int)
    CI_int = CI_int_default;
end

% Just use first kset given (best or not) as line
y_best = all_kset_tc(1,:);

% Allow multiple CI intervals - order from largest to smallest
% and then plot from lightest to darkest
assert(all(CI_int <= 1) && all(CI_int > 0),...
    "CI interval values must be between 0 and 1")
sorted_CI_ints = sort(CI_int,'descend');

best_k_color = base_color;

int_colors = repmat(base_color, length(CI_int), 1);
% Get shades multipleirs ([0.95, 0.9, 0.85, ...])

% shades_multiplier = 1 - 0.2 .* [1 : length(CI_int)]';
% shades_all = linspace(0, 1, length(CI_int) + 2);
% shades_multiplier = shades_all(end-1:-1:2)';
max_whiten_fraction = 0.3;
n_shades = length(CI_int);
shade_gap = max_whiten_fraction / (n_shades - 1);
shades_multiplier = transpose([n_shades - 1 : -1 : 0]);
if n_shades > 1
    shades_multiplier = shades_multiplier .* shade_gap;
end
% 
% int_colors = int_colors .* shades_multiplier;

int_colors = int_colors + (1 - int_colors) .* shades_multiplier;


% Plot each confidence interval
for int = 1:length(sorted_CI_ints)
    CI = CI_matrix(all_kset_tc, sorted_CI_ints(int));
    y_upper_CI = CI(2,:);
    y_lower_CI = CI(1,:);
    % Don't let lower go below 0
    y_lower_CI = max([zeros(size(y_lower_CI)); y_lower_CI]);
    CI_color = int_colors(int, :);

    CI_legend = char(strcat(num2str(sorted_CI_ints(int)*100),...
            "% CI"));

    if ~isempty(condition_name)
        CI_legend = sprintf('%s - %s',...
            condition_name,...
            CI_legend);
    end

    h_CI = fill([tc_x_vals(:)', fliplr(tc_x_vals(:)')],...
        [y_upper_CI, fliplr(y_lower_CI)],...
        CI_color,...
        'EdgeAlpha',0,...
        'DisplayName', CI_legend,...
        'FaceAlpha', 0.3);
    
    hold on

end

if plot_best_kset
    % Plot line for best kset
    kset_legend = char(sprintf("k-set %i", kset_inds(1)) );
    h_best = plot(tc_x_vals, y_best,...
        'Color', best_k_color,...
        'LineWidth', 2,...
        'DisplayName', kset_legend);
else
    h_best = NaN;
end

end