function [x3,P3]=ksmooth(x1,x2,x3,P1,P2,P3,f)
%[x3,P3]=ksmooth(x1,x2,x3,P1,P2,P3,f)
%Kalman Smoother
%given specified matrices, ksmooth computes smoothed estimates of state
%variable and its covariance matrix

 as=P1*f'/P2;
 P3=P1+as*(P3-P2)*as';
 x3=x1+as*(x3-x2); %smooth = update + weight * (smoothed from t-1 - f'* prediction)
 %in original file there is: x3=x1+as*(x3-x2) but that yields the same
 %results since x2=f*x1+z


end

 