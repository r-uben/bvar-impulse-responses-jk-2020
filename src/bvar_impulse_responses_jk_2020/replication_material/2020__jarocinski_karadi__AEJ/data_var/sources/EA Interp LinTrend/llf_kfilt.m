function [x1,x2,P1,P2,v,F,K,llf]=llf(y,x1,P1,H,T,V,W,Z)
%Given starting matrices and data, this function computed predictions,
%updates and constructs log likelihood function
    
%Log likelihood function
    %@ Kalman Filter Procedure  -- Hamilton Notation
    %   y(t) = H'x(t)+e(t)
    %   x(t) = T x(t-1) + z(t)+ v(t)
    % 
    %   var(v(t)) = W;
%P1 refers to P_t-1 (i.e. covarianc ematrix at t-1), P2=P_t|T-1 (i.e. prediction of covarianc ematrix)
%F is covariance matrix of the prediciton of y_t

%prediction of x:
x2=T*x1+Z; 
P2=T*P1*T'+W; %covar of x2
%prediction error (compared to actual observation)
v=y-H'*x2;

%prediction of y:
y2=H'*x2;
F=H'*P2*H+V; %covar of y
Finv=F\eye(size(F,1));

%Kalman gain
%k=P2*h*(inv(F));
k=P2*H/F;
K=T*P2*H/F;
%update of x1
x1=x2+k*v;
%update of p1
P1=(eye(size(P1,1))-k*H')*P2;
%@ Correct any small round-off problems @
P1=0.5*(P1+P1');

%log-likelihood function 
llf=+0.5*log(det(F))+0.5*v'/F*v; %signs are reversed as we supply negative log-lik

end

