load('SimulationResults/GeneralTrajectoryLearning/04_16_2017__08_44_21/cost_data_for_simnumber_1000_curvature_dx1000_24.mat')
load('SimulationResults/GeneralTrajectoryLearning/04_16_2017__08_44_21/std_cost_data_for_simnumber_1000_curvature_dx1000_24.mat')

cost = Inf;
for i = 1:243
    disp(i)
    a = standardized_cost_matrix(:,i);
    for j = 1:243
        b = standardized_cost_matrix(:,j);
        for k = 1:243
            c = standardized_cost_matrix(:,k);
            local_cost = min(a,min(b,c));
            if (sum(local_cost) < cost)
               cost = sum(local_cost);
               i1 =i;
               i2 =j;
               i3 =k;
            end
        end
    end
end