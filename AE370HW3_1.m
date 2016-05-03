%% AE 370 HW 3
%  Max Feinberg
%% HW3
function hw3
close all; clear all; clear figure;
 t_0 = 0.0;
 t_f = 8.0;
 y_0 = [0.2, 0.0, 0.0].';
 h = 0.01;
 t = [t_0 : h : t_f];
 y_fe = forward_euler(t, y_0, h, @yprime);
 y_be = backward_euler(t, y_0, h, @yprime);
 y_fe = y_fe(1, :);
 y_be = y_be(1, :);
 y_a = 1.02 * exp(-4.0 * t) - 2.76667 * exp(-3.0 * t) ...
 + 2.00769 * exp(-2.0 * t) - 0.0610256 * cos(3.0 * t) ...
 - 0.0682051 * sin(3.0 * t);
 figure(1)
 plot(t, y_a, 'r', t, y_fe, 'g', t, y_be, 'b', 'LineWidth', 2)
 set(gca, 'FontSize', 18)
 title(sprintf('Solutions using h = %4.2f ', h))
 xlabel('t')
 ylabel('y(t)')
 legend('Analytical', 'Forward Euler', 'Backward Euler')
 print(gcf, '-depsc2', 'hw5p1-solns')
end
%% forward euler
function y = forward_euler(t, y_0, h, yp)
 %% Integrate the ODE given by yp using the forward Euler method.
 %% Fill in required function

L = size(t);
y = zeros(1,L(1,2)-1);
w = zeros(1,L(1,2)-1);
z = zeros(1,L(1,2)-1);

y(1) = y_0(1);
w(1) = y_0(2);
z(1) = y_0(3);

for k = 1 : L(1,2)-1

    temp = yp(t(k), [y(k), w(k), z(k)]);
    z(k+1) = z(k) + h * temp(3);
    w(k+1) = w(k) + h * z(k);
    y(k+1) = y(k) + h * w(k);

end

end
%% backward_euler
function y = backward_euler(t, y_0, h, yp)
 %% Integrate the ODE given by yp using the backward Euler method.
 %% Fill in required function
tol = 10^-6;
L = size(t);
y = zeros(1,L(1,2));
y_c = 0;%zeros(1,L(1,2));
y_p = 0;%zeros(1,L(1,2));
w = zeros(1,L(1,2));
z = zeros(1,L(1,2));

y_c = y_0(1);
y(1) = y_0(1);
y_p = y_0(1);

w(1) = y_0(2);
z(1) = y_0(3);
for k = 1 : L(1,2)-1
    temp = yp(t(k), [y(k), w(k), z(k)]);
    z(k+1) = z(k) + h * temp(3);
    w(k+1) = w(k) + h * z(k);
    y_c = y(:,k) + h * w(k);
    
    y_p = y(:,k);
    countMax = -1;
    count = 0;
    while norm(y_c-y_p) > tol
       y_p = y_c;
       temp2 = yp(t(k+1), [y_p, w(k), z(k)]);
       z(k+1) = z(k) + h * temp2(3);
       w(k+1) = w(k) + h * z(k+1);
       y_c = y(k) + h * w(k+1);
       count = count + 1;
    end
    if count > countMax
       countMax = count; 
    end
    fprintf('CountMax: %f\n', countMax);
    y(:, k+1) = y_c;
end

end
%% yprime
function yp = yprime(t, y)
 %% Evaluate the derivative y'(t, y).
 yp = zeros(3, 1);
 yp(1) = y(2);
 yp(2) = y(3);
 yp(3) = - 24.0 * y(1) - 26.0 * y(2) - 9.0 * y(3) + 7.0 * sin(3.0 * t);
end 