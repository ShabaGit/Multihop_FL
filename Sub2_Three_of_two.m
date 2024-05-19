yalmip('clear')
clc;

s_l = [];

sum_m = 0; %latency for the leaf node
Lm = 20;
Cm = 1000;
Dm = 5000;
bm = 1.67e6; %20.0e6;
gm = 10^(-10/10);
n0 = 10^(-174/10);
s = 26*200e5;
pm = 10^(40/10);
fm = 2.0e9; %1.0e9;
sum_m = sum_m + ((Lm*Cm*Dm)/fm)+((s)/(bm*log2(1+((pm*gm)/(bm*n0)))));
disp(sum_m)


iter = 7; %number of SCA iteration
N = 4; %number of relay node in the route
c = [1000, 1000, 1000, 1000]; %, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000];
d = [5000, 5000, 5000, 5000]; %, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000];
g = [10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10)]; %, 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10), 10^(-10/10)];
for n = 1:1:N
    latency_values = [];
    fn_values = [];
    pn_values = [];
    yn_values = [];
    fn = sdpvar(1); %define control variable
    pn = sdpvar(1);
    yn = sdpvar(1);
    p_n_i = 1;
    y_n_i = 15;
    Pn = 10^(10/10);
    Fn = 2.0e9;
    Ln = 15;
    Cn = c(n);
    Dn = d(n);
    zeta_n = 1.0e-28;
    E_n_max = 30;
    s = 2*400e5;
    bn = 1.67e6; %15.0e6;
    n0 = 10^(-174/10);
    gn = g(n);
    for i = 1:1:iter
        cons= []; %array for storing constraints
        cons = [cons, fn >= 0];
        cons = [cons, fn <= Fn];
        cons = [cons, pn >= 0];
        cons = [cons, pn <= Pn];
        expr1 = ((n+1) * s * log(2)) / (bn);
        expr2 = yn*(log(1 + ((p_n_i * gn) / (bn * n0))) + ((p_n_i * gn) / ((pn * gn) + (bn * n0))) - ((((p_n_i * gn)^2) / ((p_n_i * gn) + (bn * n0))) * (1 / (pn * gn))));
        expr3 = Ln * zeta_n * Cn * Dn * (fn^2);
        expr4 = (1/2) * (p_n_i / y_n_i) * (yn^2);
        expr5 = (1/2) * (y_n_i / p_n_i) * (pn^2);
        lhs = expr3 + expr4 + expr5;
        cons = [cons, lhs <= E_n_max];
        cons = [cons, expr1 <= expr2];
        %disp(cons)
        expr6 = ((Ln*Cn*Dn)/(fn)); 
        Objective = expr6 + yn + sum_m;
        options = sdpsettings('verbose',1,'debug',1,'showprogress',1);
        sol= solvesdp(cons,Objective,options);
        solution = value(Objective);
        sol.info
        fprintf('fn is: %.4f\n', value(fn))
        fprintf('pn is: %.4f\n', value(pn))
        fprintf('yn is: %.4f\n', value(yn))
        fprintf('latency is: %.4f\n', value(Objective))
        %fprintf('constraint size is: %.4f\n', length(Constraints))
        latency_values(i) = value(Objective);
        fn_values(i) = value(fn);
        pn_values(i) = value(pn);
        yn_values(i) = value(yn);
        p_n_i = value(pn);
        y_n_i = value(yn);
        %fprintf('p_m_i is: %.4f\n', value(p_m_i))
        s_l = [s_l, (value(Objective)-sum_m)];
    end
    figure;
    i_values = 1:iter;
    subplot(3, 1, 1);
    plot(i_values, latency_values);
    ylabel('Latency');
    hold on;
    subplot(3, 1, 2);
    plot(i_values, fn_values);
    ylabel('fn');
    hold on;
    subplot(3, 1, 3);
    plot(i_values, pn_values);
    ylabel('pn');
    %disp(yn_values)
end


disp("Here are the s_l values:")
disp(s_l)

% Splitting into smaller arrays
smaller_arrays = reshape(s_l, 7, []);

% Result array to store the final sums
results = zeros(1, 7);

% Iterate over each position in the smaller arrays
for i = 1:7
    % Initialize the sum for this position
    position_sum = 0;
    
    % Iterate over each smaller array
    for j = 1:size(smaller_arrays, 2)
        % Add the element at the current position to the sum
        position_sum = position_sum + smaller_arrays(i, j);
    end
    
    % Add b to the sum
    position_sum = position_sum + sum_m;
    
    % Assign the result to the results array
    results(i) = position_sum;
end

% Print the results
disp("Here are the system latency values:")
disp(results);

% Plot results versus q
figure;
q = 1:7;
plot(q, results)
xlabel('Iteration')
ylabel('Latency (s)')
title('Latency versus Iteration')
grid on
