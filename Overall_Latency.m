yalmip('clear')

pn_final = [];
fn_final = [];

pm_list = [];
fm_list = [];
o2_values = [];
p_m_i_input = 0.00001;
x_m_i_input = 14;
p_n_i_input = [1, 1, 1, 1]; %, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
y_n_i_input = [15, 15, 15, 15]; %, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15];
pn_input = [10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10)]; %, 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10), 10^(10/10)];
fn_input = [0.5e9, 0.5e9, 0.5e9, 0.5e9]; %, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9, 0.5e9];

iteration = 7;

for j = 1:1:iteration
    %o1 = 0;
    o1_input = 0;

    [pm, fm, xm, o1] = SubProblem1_Function(pn_input, fn_input, p_m_i_input, x_m_i_input);
    %disp(value(pm));
    %disp(value(fm));
    %disp(value(xm));

    %o1 = o1 + ((Lm*Cm*Dm)/fm)+((s)/(bm*log2(1+((pm*gm)/(bm*n0)))));

    o1_input = o1_input + o1;
    pm_input = value(pm);
    fm_input = value(fm);
    p_m_i_input = value(pm);
    x_m_i_input = value(xm);
    pm_list(j) = value(pm);
    fm_list(j) = value(fm);
    
    %disp("ggggggggg")
    %disp(y_n_i_input) % Check the size of y_n_i_input
    %disp(y_n_i_input) % Display the values of y_n_i_input
    
    [pn_array, fn_array, yn_array, o2] = SubProblem2_Function(pm_input, fm_input, p_n_i_input, y_n_i_input, o1_input);
    disp(pn_array);
    disp(fn_array);
    disp(yn_array);
    
    o2_values(j) = value(o2);
    pn_input = pn_array;
    fn_input = fn_array;
    p_n_i_input = pn_array;
    y_n_i_input = yn_array;
    pn_final = [pn_final, pn_array];
    fn_final = [fn_final, fn_array];
end
fprintf('Here are pn values.\n');
disp(pn_final)
fprintf('Here are fn values.\n');
disp(fn_final)

fprintf('Here are pm values.\n');
disp(pm_list)
fprintf('Here are fm values.\n');
disp(fm_list)
fprintf('Here are o2 values.\n');
disp(o2_values)
%j_values = 1:iteration;
%plot(j_values, o2_values);
%ylabel('Latency');