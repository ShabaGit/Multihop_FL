function [pm, fm, xm, o1] = SubProblem1_Function(pn_input, fn_input, p_m_i_input, x_m_i_input)

%yalmip('clear')
clc;

o1 = 0;

sum = [];

N1 = 2; %latency for 2 relay nodes (route 1)
sum1 = 0;
g = [10^(-30/10), 10^(-20/10)]; %, 10^(-10/10)];
c = [10000, 10000]; %, 10000];
d = [500, 500]; %, 500];
for n = 1:1:N1
    s = 100e3;
    n0 = 10^(-174/10);
    bn = 1.67e6; %5.0e6; 
    Ln = 10;
    Cn = c(n);
    Dn = d(n);
    fn = fn_input(n);
    pn = pn_input(n);
    gn = g(n);
    sum1 = sum1 + ((Ln*Cn*Dn)/fn)+(((n+1)*s)/(bn*log2(1+((pn*gn)/(bn*n0)))));
end
%fprintf('relay sum for route 1 is: %.4f\n', sum1);
sum = [sum, sum1];

N2 = 3; %latency for 3 relay nodes (route 2)
sum2 = 0;
g = [10^(-30/10), 10^(-20/10), 10^(-10/10)]; %, 10^(-30/10), 10^(-30/10), 10^(-20/10), 10^(-10/10), 10^(-30/10), 10^(-30/10), 10^(-20/10), 10^(-10/10), 10^(-30/10), 10^(-30/10), 10^(-20/10), 10^(-10/10), 10^(-30/10), 10^(-30/10), 10^(-20/10), 10^(-10/10), 10^(-30/10)];
c = [10000, 10000, 10000]; %, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000];
d = [500, 500, 500]; %, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500];
for n = 1:1:N2
    s = 100e3;
    n0 = 10^(-174/10);
    bn = 1.67e6; %5.0e6; 
    Ln = 30;
    Cn = c(n);
    Dn = d(n);
    fn = fn_input(n);
    pn = pn_input(n);
    gn = g(n);
    sum2 = sum2 + ((Ln*Cn*Dn)/fn)+(((n+1)*s)/(bn*log2(1+((pn*gn)/(bn*n0)))));
end
%fprintf('relay sum for route 2 is: %.4f\n', sum2);
sum = [sum, sum2];

N3 = 4; %latency for 4 relay nodes (route 3)
sum3 = 0;
g = [10^(-30/10), 10^(-20/10), 10^(-10/10), 10^(-30/10)];
c = [10000, 10000, 10000, 10000];
d = [500, 500, 500, 500];
for n = 1:1:N3
    s = 100e3;
    n0 = 10^(-174/10);
    bn = 1.67e6; %5.0e6; 
    Ln = 10;
    Cn = c(n);
    Dn = d(n);
    fn = fn_input(n);
    pn = pn_input(n);
    gn = g(n);
    sum3 = sum3 + ((Ln*Cn*Dn)/fn)+(((n+1)*s)/(bn*log2(1+((pn*gn)/(bn*n0)))));
end
%fprintf('relay sum for route 3 is: %.4f\n', sum3);
sum = [sum, sum3];
%disp(sum);
sum_max = max(sum);
fprintf('maximum sum is: %.4f\n', sum_max);
sum_max_index = find(sum == max(sum));
fprintf('maximum sum index is: %.4f\n', sum_max_index);


Lm_a = [5, 5, 5];
Cm_a = [1000, 1000, 1000];
Dm_a = [50000, 50000, 50000];
bm_a = [1.67e6, 1.67e6, 1.67e6]; %[20.0e6, 20.0e6, 20.0e6];
gm_a = [10^(-25/10), 10^(-25/10), 10^(-25/10)];

j = sum_max_index;
p_m_i = p_m_i_input; 
x_m_i = x_m_i_input;
Fm = 2e9; %3e9;
Pm = 10^(12/10);
E_m_max = 100;
zeta_m = 1.0e-28;
Lm = Lm_a(j);
Cm = Cm_a(j);
Dm = Dm_a(j);
s = 2*500e6; %100e3;
bm = bm_a(j);
gm = gm_a(j);
n0 = 10^(-174/10);


iter = 1; %number of SCA iteration

latency_values = [];
fm_values = [];
pm_values = [];
xm_values = [];

fm = sdpvar(1); %define control variable
pm = sdpvar(1);
xm = sdpvar(1);


for i = 1:1:iter
    cons= []; %array for storing constraints
    cons = [cons, fm >= 0];
    cons = [cons, fm <= Fm];
    cons = [cons, pm >= 0];
    cons = [cons, pm <= Pm];
    %cons = [cons, xm >= 0];
    expr1 = (s * log(2)) / (bm);
    expr2 = xm*(log(1 + ((p_m_i * gm) / (bm * n0))) + ((p_m_i * gm) / ((pm * gm) + (bm * n0))) - ((((p_m_i * gm)^2) / ((p_m_i * gm) + (bm * n0))) * (1 / (pm * gm))));
    expr3 = Lm * zeta_m * Cm * Dm * (fm^2);
    expr4 = (1/2) * (p_m_i / x_m_i) * (xm^2);
    expr5 = (1/2) * (x_m_i / p_m_i) * (pm^2);
    lhs = expr3 + expr4 + expr5;
    cons = [cons, lhs <= E_m_max];
    cons = [cons, expr1 <= expr2];
    %cons = [cons, ((Lm*zeta_m*Cm*Dm*(fm^2))+((1/2)*(p_m_i/x_m_i)*(xm^2))+((1/2)*(x_m_i/p_m_i)*(pm^2))) <= E_m_max];
    %cons = [cons, ((s*log(2))/(bm*xm)) <= (((log(1+((p_m_i*gm)/(bm*n0))))+((p_m_i*gm)/((pm*gm)+(bm*n0)))-(((p_m_i*gm)^2/((p_m_i*gm)+(bm*n0)))*(1/pm*gm))))];
    %disp(cons)
    expr6 = ((Lm*Cm*Dm)/(fm)); 
    Objective = expr6 + xm + sum(j);
    options = sdpsettings('verbose',1,'debug',1,'showprogress',1);
    sol= solvesdp(cons,Objective,options);
    solution = value(Objective);
    sol.info
    
    %fprintf('fm is: %.4f\n', value(fm))
    %fprintf('pm is: %.4f\n', value(pm))
    %fprintf('xm is: %.4f\n', value(xm))
    %fprintf('latency is: %.4f\n', value(Objective))
    %fprintf('constraint size is: %.4f\n', length(Constraints))
    latency_values(i) = value(Objective);
    fm_values(i) = value(fm);
    pm_values(i) = value(pm);
    xm_values(i) = value(xm);
    %p_m_i = value(pm);
    %x_m_i = value(xm);
    %fprintf('p_m_i is: %.4f\n', value(p_m_i))
end

o1 = o1 + ((Lm*Cm*Dm)/fm)+((s)/(bm*log2(1+((pm*gm)/(bm*n0)))));

%{
i_values = 1:iter;
subplot(3, 1, 1);
plot(i_values, latency_values);
ylabel('Latency');
hold on;
subplot(3, 1, 2);
plot(i_values, fm_values);
ylabel('fm');
hold on;
subplot(3, 1, 3);
plot(i_values, pm_values);
ylabel('pm');

disp(xm_values)
%}
end

