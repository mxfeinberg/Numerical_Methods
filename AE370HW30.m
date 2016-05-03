%%
% AE 370 HW #3
% @author Max Feinberg
% @date 2/24/16
% @version 1.0
% Root finding
%% HW3
function hw3
close all; clear all;
 tic
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
 hold on;
 plot(t, y_a, 'r', t, y_fe, 'g', t, y_be, 'b', 'LineWidth', 2)
 set(gca, 'FontSize', 18)
 title(sprintf('Solutions using h = %4.2f ', h))
 xlabel('t')
 ylabel('y(t)')
 legend('Analytical', 'Forward Euler', 'Backward Euler')
 print(gcf, '-depsc2', 'hw5p1-solns')
 fprintf(['Both the Forward and Backward Euler schemes deviate from the analytical solution\n'...
         'at the local extrema of the function. The forward Euler scheme tends to overestimate\n' ...
         'in terms of absolute value relative the analytical solution while the backward Euler scheme\n'...
         'tends to underestimate in terms of absolute value relative to the analytical solution.\n']);
 toc
end
%% forward euler
function y = forward_euler(t, y_0, h, yp)
 % Integrate the ODE given by yp using the forward Euler method.
 % Fill in required function
L = length(t); % Used to determine looping bounds
y = zeros(3,L); %hold all y values (includes derivatives)
y(:,1) = y_0; %initial condition
for k = 1 : L-1
    y(:,k+1) = y(:,k) + h * yp(t(k), y(:,k));
end

end
%% backward_euler
function y = backward_euler(t, y_0, h, yp)
% Integrate the ODE given by yp using the backward Euler method.
tol = 10^-6; % tolerance value
L = length(t);
y = zeros(3,L);
y(:,1) = y_0; 
for k = 1 : L -1
    y_c  = y(:,k) + h *yp(t(k), y(:,k));
    y_p = y(:,k); %utilizes forward euler scheme for the first iteration
    while norm(y_c-y_p) > tol % in place iteration
       y_p = y_c;
       y_c = y(:,k) + h * yp(t(k+1), y_p);
    end
    y(:, k+1) = y_c;
end
end
%% yprime
function yp = yprime(t, y)
 % Evaluate the derivative y'(t, y).
 yp = zeros(3, 1);
 yp(1) = y(2);
 yp(2) = y(3);
 yp(3) = - 24.0 * y(1) - 26.0 * y(2) - 9.0 * y(3) + 7.0 * sin(3.0 * t);
end 