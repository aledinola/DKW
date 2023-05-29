%% Smooth IRFs

y = irf_nogrant.exit;

T = length(irf_nogrant.exit);

x = (1:T)';

n_deg = 3;
coef = polyfit(x(2:T),y(2:T),n_deg);

y_fit1 = polyval(coef,x(2:T));

y_fit = [y(1);y_fit1];

figure
plot(x,y,x,y_fit)
legend('raw','fitted')

