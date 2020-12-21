j = 0:5;
x = pi*j/10;
t1=0;
ind = 1;
for i = ind:6
    t1 = t1 + x(i)^2;
end

t2=0;
for i = ind:6
    t2 = t2 + x(i);
end

t3=0;
for i = ind:6
    t3 = t3 + x(i)*(1-cos(x(i)));
end

t4=0;
for i = ind:6
    t4 = t4 + (1-cos(x(i)));
end
%%
syms a b 
eqn1 = t1*a - t2*b - t3 == 0;
eqn2 = t2*a - (6-ind+1)*b - t4 == 0;

sol = solve([eqn1, eqn2], [a, b]);
aSol = sol.a;
bSol = sol.b;
a_val = round(aSol,4);
b_val = round(bSol,4);

%%
f=0;
for i = 1:6
    f = f + 1/12*((max(0,a_val*x(i)-b_val)-(1-cos(x(i))))^2);
end
f_val = round(f,8);

b_low = round(a_val*pi*(ind-2)/10,4);
b_high = round(a_val*pi*(ind-1)/10,4);

%%
f=0;
for i = 1:5
    f = f + 1/12*((1-cos(x(i)))^2);
end
f_val = round(f,8);

%%
f=0;
for i = 1:6
    f = f + 1/12*((1-cos(x(i)))^2);
end
f_val = round(f,8);

%% 
