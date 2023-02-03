function zz_write_elasticity_param_coeff_function(model, dv_dk_str)

% Check for and remove previous file
filename = strcat('elast_param_coeff_',model.model_name,'.m');
if exist(filename, 'file') == 2
    fclose('all');
    delete(filename)
end

if ~exist('dv_dk_str', 'var')

    [~,~,dv_dk_str] = zz_create_elasticity_coefficient_matrix(model);
end
% n_eqns = length(dv_dx);

% diary(filename)
% 
% fprintf('function dv_dx = elasticity_coeff_')
% fprintf(model.model_name)
% fprintf('(x,k,Km)')
% fprintf(1,'\n');
% fprintf(1,'\n');
% fprintf('dv_dx = zeros(length(k),length(x));')
% fprintf(1,'\n');
% fprintf(1,'\n');
%  
%  
% for i = 1:n_eqns
%     fprintf(dv_dx{i,1})
%     fprintf(1,'\n');    
% end
% 
% fprintf(1,'\n');
% fprintf('end')
% fprintf(1,'\n');
% 
% diary off

all_lines = [sprintf('function dv_dk = elast_param_coeff_%s(x,k,Km,n_params,Keq)', model.model_name);...
    '';...
    'dv_dk = zeros(length(k),n_params);';...
    '';...
    dv_dk_str;...
    '';...
    'end'];

fid = fopen(filename,'wt');
fprintf(fid, '%s\n',all_lines{:});
fclose(fid);

end