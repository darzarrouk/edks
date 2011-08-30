x = [ 1 2 4 7 ]'; %location of nearest known values
A = [ ones(4,1), x, x.^2, x.^3 ];

y = [ 5 5 7 4 ]'; %known value of function at x

c = A\y; %cubic coefficients

xn = [0:0.01:10];

yn = c(1) + c(2)*xn + c(3)*xn.^2 + c(4)*xn.^3;

plot(x,y,'o',xn,yn,'-')
