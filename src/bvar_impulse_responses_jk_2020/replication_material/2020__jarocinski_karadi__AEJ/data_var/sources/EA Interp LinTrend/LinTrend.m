function [trend]=LinTrend(X);
%Fits linear trend using least-squares method
T=[1:size(X,1)]';
XX=X(~isnan(X)); %XX and TT drop NaN values and account for them
TT=T(~isnan(X));
fun=@(x)sum((XX-(x(1)+x(2)*TT)).^2); %create function to supply to optim. function
%function is a sum of sqaured residuals, which we want to minimize
x0=[1, 1]; %starting values - just to let Matlab know what form of answer we expect
[x, fval]=fminunc(fun,x0); %find minimum of a function

trend=x(1)+x(2)*T; %create fitted series based on optimized parameters
end



