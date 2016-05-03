close all; clear all; clear graph;
%%
currValPos = 0;
currValNeg = 0;
ref = @(x) log(1+x);
test = 0.5;
test2 = -0.5;
errorPos1 = 1:1000;
errorPos2 = 1:1000;
for i=1:1000
     currValPos = currValPos + ((test^i)*((-1)^(i-1))*(1/i));
     currValNeg = currValNeg + ((test2^i)*((-1)^(i-1))*(1/i));
     errorPos1(i) = abs(((currValPos - ref(test))/ref(test)));
     errorPos2(i) = abs(((currValNeg - ref(test2))/ref(test2)));
end
figure(1)
loglog(1:1000,errorPos1,'r-',1:1000,errorPos2,'b-');
grid on
xlabel('Number of Terms')
ylabel('Relative Error')
xlim([1,1000])
ylim([10e-17,1])
%%
close all; clear all;
currValPos = 0;
ref = @(x) sin(x);
test = pi/4;
errorPos1 = 1:1000;
nValues = 1:5;
exponentCheck = -6;
numFound = 1;
testStuff = 1:2;
%test = 0:2*pi/49:2*pi;
for i=1:1000
     currValPos = currValPos + ((test^(2*i-1))*((-1)^(i-1))*(1/factorial(2*i-1)));
     testStuff(i) = currValPos;
     errorPos1(i) = abs(((currValPos - ref(test))/ref(test)));
     if(errorPos1(i) <= (10^(exponentCheck)))
         nValues(numFound) = i;
         numFound = numFound + 1;
        exponentCheck = exponentCheck - 2;
     end
     if(exponentCheck == -16)
         break
     end
       
end
figure(2)
plot(0:2*pi/50:2*pi,testStuff);
% loglog(1:1000,errorPos1,'r-');
% grid on
% xlabel('Number of Terms')
% ylabel('Relative Error')
% xlim([1,1000])
% ylim([10e-17,1])