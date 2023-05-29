function make_plots_ss(sol,distribS,b_grid,par)

mu_x = squeeze(sum(distribS.mu,[1,2]));
mu_active_x = squeeze(sum(distribS.mu_active,[1,2]));

mu_k = squeeze(sum(distribS.mu,[2,3]));
mu_active_k = squeeze(sum(distribS.mu_active,[2,3]));

mu_b = squeeze(sum(distribS.mu,3));
mu_active_b = squeeze(sum(distribS.mu_active,3));

% marginal distribution of x conditional on k
mu_x_k = squeeze(sum(distribS.mu,2));
mu_active_x_k = squeeze(sum(distribS.mu_active,2));

figure
plot(par.x_grid,mu_x)
hold on 
plot(par.x_grid,mu_active_x)
legend('mu','mu_active')
xlabel('Productivity, x')
title('Marginal distribution of x')



k_c = 2;
figure
plot(par.x_grid,squeeze(mu_x_k(k_c,:)))
hold on
plot(par.x_grid,squeeze(mu_active_x_k(k_c,:)))
legend('mu','mu_active')
xlabel('Productivity, x')
title('Marginal distribution of x, conditional on k')



figure
plot(par.k_grid,mu_k)
title('Marginal distribution of k')


figure
k1 = 1; %round(nk/2);
k2 = 20;
k3 = 40;
k4 = par.nk;
plot(squeeze(b_grid(k1,:)),squeeze(mu_b(k1,:)/sum(mu_b(k1,:),"all")))
hold on
plot(squeeze(b_grid(k2,:)),squeeze(mu_b(k2,:)/sum(mu_b(k2,:),"all")))
hold on
plot(squeeze(b_grid(k3,:)),squeeze(mu_b(k3,:)/sum(mu_b(k3,:),"all")))
hold on
plot(squeeze(b_grid(k4,:)),squeeze(mu_b(k4,:)/sum(mu_b(k4,:),"all")))
legend('k low','k medium low','k medium high','k high')
xlabel('Debt, b')
title('Marginal distribution of b, conditional on k')

end

