function [Qtrend, Mtrend]=LinTrend_qm(X);
%Fits linear trend using least-squares method
T=[1:size(X,1)]';
Tm=[1:3*size(X,1)]'; %T used for monthly trend
fun=@(x)sum((X-(x(1)+x(2)*T)).^2); %create function to supply to optim. function
%function is a sum of sqaured residuals, which we want to minimize
x0=[1, 1]; %starting values - just to let Matlab know what form of answer we expect
[x, fval]=fminunc(fun,x0); %find minimum of a function

Qtrend=x(1)+x(2)*T; %create fitted series based on optimized parameters

T_m=[1:3*size(X,1)]';
Mtrend=x(1)+x(2)/3*Tm; %create monthly trend
end
