close all; clear all; clear figure;
%% Problem 1
%AE 370 HW #1
% @author Max Feinberg
% @date 2/10/16
% @version 1.0
% Numerical Integration and Differentation
%%Problem 1
f = @(x) sinh(x);
Problem1_Check = integral(f,2,6);
%Rectangle
a = 2; b = 6;
h = (b-a)/3;
Rectangle_Rule_Approximation = h * (f(8/3)+f(12/3)+f(16/3));
Rectangle_Error = (Rectangle_Rule_Approximation -  Problem1_Check) / Problem1_Check * 100;
%Trapezoid
Trapezoidal_Rule_Approximation = (2/2)*(f(2)+2*f(4)+f(6));
Trap_Error = (Trapezoidal_Rule_Approximation -  Problem1_Check) / Problem1_Check * 100;
%Simpson
Simpson_Rule_Approximation = (h/2)*(f(2)+4*f(4)+f(6));
Simpson_Error = (Simpson_Rule_Approximation -  Problem1_Check) / Problem1_Check * 100;
%Gauss Quadrature
zeta1 = 0; w_1 = 8/9;
zeta2 = sqrt(3/5); w_2=5/9;
zeta3 = -1*sqrt(3/5);
GQ_Approximation = (b-a)*.5*(f(((b-a)*zeta1 + (a+b))/2)*w_1 + f(((b-a)*zeta2 + (a+b))/2)*w_2 + f(((b-a)*zeta3 + (a+b))/2)*w_2);
GQ_Error = (GQ_Approximation -  Problem1_Check) / Problem1_Check * 100;
fprintf(['Problem 1\n'...
        '\nExact Integral Value: %f \n'...
        '\nRectangle Rule Value: %f Percentage Error: %f%% \n'...
        '\nTrapezoidal Value: %f Percentage Error: %f%% \n'...
        '\nSimpson''s Rule Value: %f Percentage Error: %f%% \n'...
        '\nGauss Quadrature: %f Percentage Error: %f%% \n'], Problem1_Check, Rectangle_Rule_Approximation, Rectangle_Error,...
        Trapezoidal_Rule_Approximation, Trap_Error, Simpson_Rule_Approximation, Simpson_Error, GQ_Approximation, GQ_Error);
fprintf(['With only 3 function evaluations, the Gauss Quadrature method of numerical integration is clearly the most accurate\n'...
         'with a percentage error of only -0.126%%.  As is expected, the Rectangle Rule Approximation is more than twice as accurate\n'...
         'as the Trapezoidal Rule Approximation.  The Simpson''s rule approximation also yields a more accurate result than the\n' ...
         'Rectangle Rule and Trapezoidal Rule approximations.\n\n']);
     
%% Problem 2
fprintf('Problem 2\n');
f2 = @(x) cos((x.^3)+3.*x)./exp(x.^2);
Problem2_Check = integral(f2,0, 2*pi);
TrapezoidalValues = zeros(30,1);
k = 0;
MAX_INTERVAL = 30;
RectangleValues = zeros(2,MAX_INTERVAL);
TrapValues = zeros(2,MAX_INTERVAL);
while k < MAX_INTERVAL
    k = k+1;
    w = (2*pi)/(k);
    Ir = 0;
    It = 0;
    j = 0;
    while j < k
        j = j + 1;
        x = (j-1)*w+(w/2);
        z0 = x - (w/2);
        z1 = x + (w/2);
        Ir = Ir + f2(x)*w;
        It = It + (f2(z0) + f2(z1))*.5*w; 
    end
    RectangleValues(1,k) = k;
    RectangleValues(2,k) = Ir;
    TrapValues(1,k) = k;
    TrapValues(2,k) = It;
end
figure(1)
plot(RectangleValues(1,:),RectangleValues(2,:),'b.--',TrapValues(1,:),TrapValues(2,:),'r.-');
title('Numerical Integral Value vs. Number of Intervals');
legend('Rectangle', 'Trapezoidal');
ylabel('Integral value');
xlabel('Number of Intervals');
grid on;
fprintf(['Both the Rectangle Rule and Trapezoidal Rule approximations begin to converge\n' ...
         'at approximately 15 intervals. At low interval values, the Trapezoidal Rule tends\n'...
         'to overestimate the integral for this function while the Rectangle rule tends to \n'...
         'underestimate the integral.\n\n']);
%% Problem 3
fprintf('Problem 3\n');
f3 = @(x,y) x.*sin(x.^2)+log(2 + y);
Problem3_Check = integral2(f3, -1, 1, -1, 1);
% 1 x 1
GQ1 = 4 * f3(0,0);
GQ1_Error = (GQ1 -  Problem3_Check) / Problem3_Check * 100;
% 2 x2 
eta = sqrt(3)^-1;
zeta = sqrt(3)^-1;
GQ2 = 1 * (f3(eta, zeta)) + 1 * (f3(-1*eta, zeta)) + 1 * (f3(eta, -1*zeta)) + 1 * (f3(-1*eta, -1 * zeta));
GQ2_Error = (GQ2 -  Problem3_Check) / Problem3_Check * 100;
%3 x 3
w1 = 8/9;
w2 = 5/9;
eta = sqrt(3/5);
zeta = sqrt(3/5);
GQ3 = w1 * w1 * f3(0,0) + w1 * w2*f3(eta, 0) + w1 * w2 * f3(-1 * eta, 0) + w1 * w2 * f3(0, eta) + ...
      w1 * w2 * f3(0, -1 * eta) + w2 * w2 * f3(eta, zeta) + w2 * w2 * f3(-1 * eta, zeta) +...
      w2 * w2 * f3(eta, -1 * zeta) + w2 * w2 * f3(-1 * eta, -1 *zeta);
GQ3_Error = (GQ3 -  Problem3_Check) / Problem3_Check * 100;
fprintf(['Exact Integral Value: %f\n'...
         '\n1 * 1 Gauss Quadrature Rule Approximation: %f Percentage Error: %f%%\n' ...
         '\n2 * 2 Gauss Quadrature Rule Approximation: %f Percentage Error: %f%%\n' ...
         '\n3 * 3 Gauss Quadrature Rule Approximation: %f Percentage Error: %f%%\n\n' ...
        ], Problem3_Check, GQ1, GQ1_Error, GQ2, GQ2_Error, GQ3, GQ3_Error );
fprintf(['The Guass Quadrature approximations of the integral quickly converge to the\n'...
         'exact value as the number of sampling points increases.\n\n']);
%% Problem 4
fprintf('Problem 4\n');
f4 = @(x) exp(1 + 3.*x);
syms t;
Problem4_Check = diff(exp(1+3*t), t);
t = 2;
Problem4_Check = eval(Problem4_Check);
h1 = 0.05; h2 = 0.1;
x = 2;
CentralDifferenceApprox1 = (f4(x + h1) - f4(x-h1))/(2*h1);
CD1_Error = (CentralDifferenceApprox1 -  Problem4_Check) / Problem4_Check * 100;
CentralDifferenceApprox2 = (f4(x + h2) - f4(x - h2))/(2*h2);
CD2_Error = (CentralDifferenceApprox2 -  Problem4_Check) / Problem4_Check * 100;
a0 =(CentralDifferenceApprox2 - (2^2)*CentralDifferenceApprox1)/(1-(2^2));
R_Error = (a0 -  Problem4_Check) / Problem4_Check * 100;
fprintf(['Exact Derivative Value: %f\n'...
         '\nh = 0.05 Central Difference Approximation: %f Percentage Error: %f%%\n' ...
         '\nh = 0.1 Central Difference Approximation: %f Percentage Error: %f%%\n' ...
         '\nRichardson''s Extrapolation Scheme Approximation: %f Percentage Error: %f%%\n\n' ...
        ], Problem4_Check, CentralDifferenceApprox1, CD1_Error, CentralDifferenceApprox2, ...
          CD2_Error, a0, R_Error );
fprintf(['Beginning with a small h value of 0.05, our central difference approximation has\n'...
         'a relatively small error of 0.375%%.  However, by implementing the Richardson Extrapolation method\n'...
         'with a q value of 2, we are able to reduce our error drastically to -0.00169%%, a value that is correct\n' ...
         'to 5 significant figures.']);