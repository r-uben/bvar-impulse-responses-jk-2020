function [x1,x2,P1,P2,e,F,llf]=llf(y,x1,P1,h,f,r,q,z)

%Given starting matrices and data, this function computed predictions,
%updates and constructs log likelihood function
    
%Log likelihood function

    %function [x1,p1,x2,p2,e,ht,llf]=kfilt_llf(y,x1,p1,h,f,r,q,z)
    %proc (7)=kfilt_llf(y,x1,p1,h,f,r,q,z);

    %@ Kalman Filter Procedure  -- Hamilton Notation
    %   y(t) = h'x(t) + w(t)
    %   x(t) = f x(t-1) + z(t)+ v(t)
    % 
    %   var(w(t)) = r;
    %   var(v(t)) = q;
    % 
    %   x1 = x(t-1/t-1)  -- on input
    %   p1 = p(t-1/t-1)  -- on input

    %  Output includes innovations, variance of innovations and log-likelihood;
    %

%prediction of x:
x2=f*x1+z; 
P2=f*P1*f'+q; %covar of x2
%P2
%prediction error (compared to actual observation
e=y-h'*x2;

%prediction of y:
y2=h'*x2;
F=h'*P2*h+r; %covar of y

%Kalman gain
%k=P2*h*(inv(F));
k=P2*h/F;
%update of x1
x1=x2+k*e;
%update of p1
P1=(eye(size(P1,1))-k*h')*P2;
%@ Correct any small round-off problems @
P1=0.5*(P1+P1');

%log-likelihood function 
%llf=+0.5*log(det(F))+0.5*e'*inv(F)*e; %signs are reversed as we supply negative log-lik (but even if multiplied by (-1) the convergence does not occur)
llf=+0.5*log(det(F))+0.5*e'/F*e; %signs are reversed as we supply negative log-lik (but even if multiplied by (-1) the convergence does not occur)
%llf
end

