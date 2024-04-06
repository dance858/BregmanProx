% Running this script produces one of the figures in the paper.
clear; clc;

% Load problem instance
p = 50; run = 2 ;
load(sprintf("../results/FBPG_p=%d_m=10_run=%d.mat", p, run))
load(sprintf("../results/cvx_p=%d_m=10_run=%d.mat", p, run))

%% Primal and dual bounds
objs = FBPG_struct.objs;
lower_bound = FBPG_struct.lower_bound;
cvx_optval = cvx_struct.cvx_optval;
figure()
plot(1:length(objs), (objs - cvx_optval)/cvx_optval);
hold on
plot(1:length(lower_bound), (lower_bound - cvx_optval)/cvx_optval)
xlabel('$k$', 'Interpreter', 'Latex')
ylabel('Optimality gap', 'Interpreter', 'Latex')
grid on;
yscale_symlog;
yline(0, '--')
legend('$(f(\mathcal{X}_k) - f^\star)/f^\star$', '$(L(\mathcal{X}_k) - f^\star)/f^\star$', 'Interpreter', 'Latex')
xlim([0, 120])
