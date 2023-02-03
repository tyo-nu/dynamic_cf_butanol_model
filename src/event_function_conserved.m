% EVENT_FUNCTION_CONSERVED Custom event function for integration of ODEs.
%
% function [value,isterminal,direction] = ...
%     event_function_conserved(t,xi,CPU_t,kinetic_param,S_f_b,Sr,T,Lo,...
%     mass_balance_ode_conserved,denom_metabs,denom_Kms,Options)

% Revision history:
%{
2019-03-29: jpm
    Added header
    New args in of denom_metabs, denom_Kms, needed to pass in to
    mass_balance_ode_conserved

%}

function [value,isterminal,direction] = ...
    event_function_conserved(t,xi,CPU_t,k,S_f_b,Sr,T,Lo,...
    mass_balance_fxn,denom_metabs,denom_Kms,Options,...
    ref_metab_concs, Keq, inst_met_inds, inst_rxn_inds, S_inst, ...
    species_cmp_vols,fixed_met_subs_rxns, fixed_met_timepoints, fixed_met_concs);

if Options.event_fxn_ignore_dxdt
    Error = 1;
else
    % Calculate maximum deviation from dx = S*v = 0
    Error = norm(mass_balance_fxn(t,xi,k, ...                         
        S_f_b,Sr,T,Lo,denom_metabs,denom_Kms,ref_metab_concs,...
        Keq, inst_met_inds, inst_rxn_inds, S_inst, species_cmp_vols,...
        fixed_met_subs_rxns, fixed_met_timepoints, fixed_met_concs));
end

% Required for event function
direction = 0;                                                          

%pass in new values from Options
% Max dxdt allowed
tol = Options.event_fxn_tolerance.*10^-6;   
% Max time in integration allowed
timeout = Options.event_fxn_timeout;

time_elapsed = toc(CPU_t);

% If integration time is greater than 300s OR error drops below tolerance, end integration; 
if time_elapsed < timeout && Error>tol                                      
    value=Error;
    isterminal = 0;
    
% Avoid event not triggering if event happens on first integration step
elseif t <= 0.001                                                        
    value=Error;
    isterminal = 0;

% Terminate integration if error drops below tolerance or integration time
% surpases limit
else
    value=0;
    isterminal = 1;                                                     
end
    
end