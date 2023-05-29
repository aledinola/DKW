function F = fun_tfp_small(mu,xprod_grid,gamma2)
% fun_tfp_small computes the x_bar (i.e. optimal effective TFP in 
% the small firm sector). See memo
% misallocation.pdf written by Leo on October 21, 2022
% Inputs:
% - mu: dim(nk,nb,nx)
% - xprod_grid: dim(nx,1)
% - gamma2: scalar

nu_x = squeeze(sum(mu,[1,2]));

F = sum(xprod_grid.^(1/(1-gamma2)).*nu_x)^(1-gamma2);

end %end function
