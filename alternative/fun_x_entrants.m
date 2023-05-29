function [x0_prob] = fun_x_entrants(x_grid,x_prob,epsx,rhox,mean_x,xi)
% This function computes the x-distribution for new entrants 
% (prod shifted by xi)

nx = length(x_grid);

var_x=epsx^2/(1-rhox^2);

x0_prob = x_prob;
mean_x_exp = exp(mean_x/(1-rhox));
for x_c = 2:nx-1
    x_val  = x_grid(x_c);
    x1_val = x_grid(x_c+1);
    x0_val = x_grid(x_c-1);
    x0_prob(x_c) = exp(-(log(x_val)-log(xi*mean_x_exp))^2/2/var_x)/x_val/var_x^0.5/(2*pi)^0.5*(x1_val-x0_val)/2;
end
x0_prob(1)  = exp(-(log(x_grid(1))-log(xi*mean_x_exp))^2/2/var_x)/x_grid(1)/var_x^0.5/(2*pi)^0.5*(x_grid(2)-x_grid(1))/2;
x0_prob(nx) = exp(-(log(x_grid(nx))-log(xi*mean_x_exp))^2/2/var_x)/x_grid(nx)/var_x^0.5/(2*pi)^0.5*(x_grid(nx)-x_grid(nx-1))/2;

x0_prob = x0_prob/sum(x0_prob);

end %end function

