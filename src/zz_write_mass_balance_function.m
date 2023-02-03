function zz_write_mass_balance_function(model, v_str)

% Check for and remove previous file
filename = strcat('mass_balance_ode_',model.model_name,'.m');
if exist(filename, 'file') == 2
    fclose('all');
    delete(filename)
end

if ~exist('v_str','var')
    [~,v_str] = zz_create_elasticity_coefficient_matrix(model);
end

% n_eqns = length(v_str);
% 
% diary(filename)
% 
% 
% fprintf('function dxdt = mass_balance_ode_')
% fprintf(model.model_name)
% fprintf('(S_f_b,x,k,Km)')
% fprintf(1,'\n');
% fprintf(1,'\n');
% fprintf('v = zeros(length(k));')
% fprintf(1,'\n');
% fprintf(1,'\n');
%  
%  
% for i = 1:n_eqns
%     fprintf(v{i,1})
%     fprintf(1,'\n');    
% end
% 
% fprintf('dxdt = S_f_b * v;')
% 
% fprintf(1,'\n');
% fprintf('end')
% fprintf(1,'\n');
% 
% diary off

all_lines = [sprintf('function [dxdt, v] = mass_balance_ode_%s(~,x,k,Km,S_f_b)', model.model_name);...
    '';...
    'v = zeros(length(k), 1);';...
    '';...
    v_str;...
    '';...
    'dxdt = S_f_b * v;';...
    '';...
    'end'];

fid = fopen(filename,'wt');
fprintf(fid, '%s\n',all_lines{:});
fclose(fid);

end