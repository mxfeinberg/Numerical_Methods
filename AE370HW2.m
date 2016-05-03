
%%
%AE 370 HW #2
% @author Max Feinberg
% @date 2/10/16
% @version 1.0
% Root finding
%% Main
function hw2p1
% hw2p1 AE 370 Spring 2016 Homework 2 Problem 1
close all; clear all; clear figure;
 x = linspace(0.0, 5.0, 101); % Plot f(x) and f'(x).
 f = fn(x);
 fp = fnp(x);
 figure (1)
 plot(x, f, 'k', 'LineWidth', 2)
 grid on
 set(gca, 'FontSize', 18)
 title('Given function f(x)')
 xlabel('x')
 ylabel('f(x)')
 print(gcf, '-depsc2', 'hw2p1-fn')

 x_low = 1.0;
 x_high = 4.0;
 x_tol = 1.0e-6;
 x_hat_bis = bisection(@fn, x_low, x_high, x_tol);
 x_0 = x_low + 0.5 * (x_high - x_low); % Starting guess
 f_tol = 1.0e-6;
 x_hat_nr = newton_raphson(@fn, @fnp, x_0, f_tol);
 x_1 = x_0 - 0.01 * (x_high - x_low); % Second starting guess
 x_hat_sec = secant(@fn, x_0, x_1, f_tol);

 fprintf('Bisection: x_hat = %10.6f\n', x_hat_bis)
 fprintf('Newton-Raphson: x_hat = %10.6f\n', x_hat_nr)
 fprintf('Secant method: x_hat = %10.6f\n', x_hat_sec)
end

%% function
function f = fn(x)
% Evaluate the given function.
 f = 2.0 * x.^3 + 5.875 * x.^2 - 8.625 * x - 24.75;
end
%% Derivative
function fp = fnp(x)
% Evaluate the derivative of the given function.
 fp = 6.0 * x.^2 + 11.75 * x - 8.625;
end
%% Bisection
function x = bisection(fBIS, x_low, x_high, x_tol)
% Isolate a root of f(x) using bisection.
% fill up this code here
    a=x_low;
    b=x_high;
    if fBIS(a)*fBIS(b) >=0
        sprintf( 'this interval does not contain a root')
        x = sqrt(-1);
    else
        tol=x_tol;
        %sprintf('a=%17.9f f(a)=%17.9f b=%17.9f f(b)=%17.9f\n',a,fn(a),b,fn(b))
        while b-a > tol
            c=(a+b)/2;
            if fBIS(a)*fBIS(c) > 0
                 a=c;
            else
                 b=c;
            end
        %sprintf('a=%17.9f f(a)=%17.9f b=%17.9f f(b)=%17.9f\n',a,fBIS(a),b,fBIS(b))
        end
    %sprintf(' the root is in the interval %17.9f ... %17.9f\n',a,b)
    end
    x = (a+b)/2;

end
%% Newton Raphson
function x = newton_raphson(f, fp, x_0, tol0)
% Isolate a root of f(x) using Newton-Raphson iteration.
% fill up this code here

x=x_0;
tol=tol0;
iter_max=20;
iter=0;
while abs(f(x)) > tol
    iter=iter+1;
    if iter > iter_max
        break
    end
    dx=-f(x)/fp(x);
    %sprintf(' iter= %d x=%15.10f f(x)=%20.15f dx=%20.15f\n',iter,x,f(x),dx)
    x=x+dx;
end
%sprintf(' the root is %20.15f : f(x) = %20.15f\n',x,f(x))

end

%% Secant
function x = secant(f, x_0, x_1, tol0)
% Isolate a root of f(x) using the secant method.
% fill up this code here
x0=x_0;
x1=x_1;
x = x1;
tol=tol0;
iter_max=20*2;
iter=0;
while abs(f(x)) > tol
    iter=iter+1;
    if iter > iter_max
        break
    end
    temp = x;
    x = x - f(x)*((x-x_0)/(f(x)-f(x_0)));
    x_0 = temp;
    %sprintf(' iter= %d x=%15.10f f(x)=%20.15f',iter,x,f(x));
end
%sprintf(' the root is %20.15f : f(x) = %20.15f\n',x,f(x))
end