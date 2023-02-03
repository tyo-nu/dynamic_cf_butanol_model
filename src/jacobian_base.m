% Just bare-bones jacobian calculator

function df_dx = jacobian_base(~, x, k, Km, S_f_b, elasticity_coeff_fxn)

e_ij = elasticity_coeff_fxn(x, k, Km);

if any(isnan(e_ij),'all')
    error("Something is wrong with elasticity coeff")
end

df_dx = S_f_b * e_ij;

df_dx = sparse(df_dx);

end