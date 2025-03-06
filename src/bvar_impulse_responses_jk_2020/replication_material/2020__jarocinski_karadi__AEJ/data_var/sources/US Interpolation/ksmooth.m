function [x3,p3]=ksmooth(x1,x2,x3,p1,p2,p3,f)
%[x3,P3]=ksmooth(x1,x2,x3,P1,P2,P3,f)
%Kalman Smoother
%given specified matrices, ksmooth computes smoothed estimates of state
%variable and its covariance matrix

 as=p1*f'/p2;
 p3=p1+as*(p3-p2)*as';
 x3=x1+as*(x3-x2); %smooth = update + weight * (smoothed from t-1 - f'* prediction)
 %in original file there is: x3=x1+as*(x3-x2) but that yields the same
 %results since x2=f*x1+z


 end
