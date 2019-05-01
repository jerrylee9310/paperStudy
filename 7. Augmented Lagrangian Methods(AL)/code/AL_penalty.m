close all; clear all; clc;

x = linspace(0,5,500);
f = (x-2).^2 + 1;
g = sqrt(x)-1;
figure();
plot(f); hold on; plot(g); hold off;

% Quadratic function 
    % yield slightly infeasible solution
    % make it more convex ( but too big weight on qua term make it ill-con)
penalty_weight = [1 2 4 10 100];
penalty = [];
for w = penalty_weight
    penalty = [penalty; f + w*(g.^2)];
end
figure(); 
for i = 1:5
    plot(penalty(i,:));
    hold on;
end
hold off;

lagrangian_lambda = [1 2 4 8];
lg = [];
for w = lagrangian_lambda
    lg = [lg; f + w*g];
end
figure();
for i = 1:4
    plot(lg(i,:));
    hold on;
end
hold off;