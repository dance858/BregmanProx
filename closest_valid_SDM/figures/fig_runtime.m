% Running this script produces one of the figures in the paper.
clear; clc;
%% Parse data
number_of_instances = 3;
all_p_plot = [20, 30, 40, 50, 60, 70, 80];
m = 10;
runtime_cvx = zeros(number_of_instances, length(all_p_plot));
runtime_FBPG = zeros(number_of_instances, length(all_p_plot));
runtime_FBPG_new = zeros(number_of_instances, length(all_p_plot));

iterations_IGA = zeros(number_of_instances, length(all_p_plot));
average_Newton_steps = zeros(number_of_instances, length(all_p_plot));
max_Newton_steps = zeros(number_of_instances, length(all_p_plot));

for k = 1:length(all_p_plot)
    p = all_p_plot(k);
    for instance = 1:number_of_instances
        load(sprintf("../results/cvx_p=%d_m=%d_run=%d.mat", all_p_plot(k), ...
                     m, instance), 'cvx_struct');
        runtime_cvx(instance, k) = cvx_struct.solve_time;
        
        load(sprintf("../results/FBPG_p=%d_m=%d_run=%d.mat", all_p_plot(k), ...
                     m, instance), 'FBPG_struct')
        runtime_FBPG(instance, k) = FBPG_struct.time;
        iterations_IGA(instance, k) = length(FBPG_struct.objs);
        average_Newton_steps(instance, k) = mean(FBPG_struct.stats_proj(1, :));
        max_Newton_steps(instance, k) = max(FBPG_struct.stats_proj(1, :));
    end
end

%% Run time figure IPM vs PG
semilogy(all_p_plot, mean(runtime_FBPG), '-o');
hold on
semilogy(all_p_plot, mean(runtime_FBPG_new), '-o');
hold on
semilogy(all_p_plot, mean(runtime_cvx), '-o');
grid on
legend('IGA', 'IPM', 'Location','northwest')
xlabel('$p$','Interpreter','Latex')
ylabel('Time (seconds)', 'Interpreter', 'Latex', 'fontsize', 12)