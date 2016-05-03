x = linspace(0,1, 100);
y0 = (1/(3*5)).*x.^2;
y1 = (31/(48*4)) .* x.^2 + -1*(13/(48*4)).*x.^3;
y2 = (1/24).*x.^4 + ((-29)/(144+48)).*x.^3 + (13/(48+16)).*x.^2;
plot(x, y0,'rs', x, y1,'bo',  x, y2, 'gd')
xlabel('x')
ylabel('EIw(x)')
title('EIw(x) vs x, Max A. Feinberg')
legend('ax^2', 'ax^2+bx^3', 'Exact Solution')