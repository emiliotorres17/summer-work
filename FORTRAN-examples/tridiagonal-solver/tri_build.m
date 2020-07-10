close all
clear all
clc

n   = 25;

a   = randn(1,n)
b   = randn(1,n+1)
c   = randn(1,n)
d   = randn(1,n+1)
lhs = zeros(n+1, n+1)

for i = 2:n
    lhs(i,i-1)  = a(i-1);
    lhs(i,i)    = b(i);
    lhs(i,i+1)  = c(i)
end

lhs(1,1)        = b(1);
lhs(1,2)        = c(1);
lhs(n+1, n+1)   = b(n+1);
lhs(n+1, n)     = a(n);

disp(lhs)


rhs = d'

x   = inv(lhs)*rhs

save('tri_example')
