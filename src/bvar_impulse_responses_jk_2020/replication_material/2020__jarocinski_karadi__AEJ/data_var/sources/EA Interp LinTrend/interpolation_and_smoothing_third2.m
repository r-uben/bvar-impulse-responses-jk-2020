function [tmp2]=interpolation_and_smoothing(b_st)
%y=interpolation_and_smoothing(b_st)
%given coefficients b-st for KF, the function computes smoothed estimates
%of state vector, using Durbin and Koopman algorithm

% The system:
%y_t=H*x_t+e_t
%x_t=f*x_t-1+z+v_t
%Var(e_t)=V=0
%Var(v_t)=W


%Note: function also uses global variables as inputs
 global Yreg; global Xreg2; global Y_m; global Y_q; global i; global Y_ext; global element;

 H = zeros(4,1);
 V=0;

 T = zeros(4,4);
 T(1,4) = 1;
 T(2:3,1:2) = eye(2);
 
 W = zeros(4,4);  %q is covariance matrix in the transition eq
 
 s=size(T,1); 
 w=size(Y_ext,1);
 %FFF, KKK, PPP, vvv, rrr are matrices used to store F,K,P,w and r
 FFF=zeros(size(H',1),size(H,2),w); 
 KKK=zeros(size(T,1),size(FFF,1),w); PPP=zeros(s,s,w); WWW=zeros(s,s,w);
 vvv=zeros(1,1,w);
 rrr=zeros(s,w); aaa=zeros(s,w);
 
 
 %% count number of monthly indicators from Xreg2 for a given quarterly series - store their number in numind
 numind=0;
 for j=1:size(Xreg2,2)
    if size(Xreg2(j).(element),2)>0
        numind=numind+1;
    end     
 end
 %% compute fitted values of monthly indices
 
 %Split b into components for each xreg (monthly index) separately
 index_size=0;
 for q=1:numind
     b(q).(element)=b_st(index_size+1 : index_size+(size(Xreg2(q).(element),2)+2)) ;
     index_size=index_size+size(Xreg2(q).(element),2)+2;
 end
 %fitted values:    
 Zreg=NaN(size(Xreg2(1).(element),1),numind);
 for q=1:numind
     Zreg(:,q)=Xreg2(q).(element)*[b(q).(element)(1:size(Xreg2(q).(element),2))]; %zreg1 = xreg1*bx1; each colum of Zreg contains fitted value from one set of regressors
 end
 
 %% 
 %@ Initial Values @
 tmp = max(abs(Yreg(:,i)));
 P1 = tmp*tmp*eye(size(T,1));
 x1 = zeros(size(T,1),1); %starting value for x (unobserved component) 
 
 %Now: fill the z matrix, i.e. the regression component in the transition
 %equation (using zreg1 and zreg2 from above)
 %zreg2 (based on xreg2 takes precedence - only if it is missing is zreg1 and
 %xreg1  are used
 for  im=1:size(Zreg,1) %i.e. while im is less or equal to the number of monthly observations
     
     iq = im/3; %covnerts im into quarters: since we start in January, Q1 this is valid
     %iq is integer for months 3 6 9 and 12, otherwise its a fraction
     z=NaN(4,1); %4x1 matrix of missing values 
     Fill=0; %indicates whether the position for obs im has yet been filled in by fitted values of indep regressors

     %Now fill in z matrix - external component in the transition equation
     q=size(Zreg,2); %determines which index of exog regressors is considered - for the next 'while' loop
     while Fill<1 
         if max(isnan(Zreg(im,q)))==0 %if Zreg(im,q) is non-NaN - fill in the position and change Fill=1 - terminate loop/(search)
             z = vertcat(Zreg(im,q),zeros(3,1)); 
             zzz(:,im)=z; %save zreg
             T(4,4) = b(q).(element)(size(Xreg2(q).(element),2)+1);
             fff(:,:,im)=T;
             W(4,4) = b(q).(element)(size(Xreg2(q).(element),2)+2)^2;
             Fill=1;
         else %if Zreg(im,q) is NaN - leave Fill=0 and continue search and 
              %reduce q by 1 (look in different column of Z i.e. try different
              %indep regressor)
             q=q-1;
         end
     end
     
     %Get matrices needed for smoothing
     if iq - floor(iq)== 0 %if iq is integer, i.e. if month=3,6,9 or 12
         y = Yreg(iq,i); %no change in y
         H(1) = Y_m(im,i)/(3*Y_q(iq,i)); %smt is a monthly trend from cubic spline
         H(2) = Y_m(im-1,i)/(3*Y_q(iq,i)); %smt and sqt re-weight to take trends into account, 1/3 come from averaging - 1 quarter has 3 months
         H(3) = Y_m(im-2,i)/(3*Y_q(iq,i));
         HHH(:,im)=H;
         
         %feed data to llf function to get likelihood function for this
         %time period to get matrices needed
         [x1,x2,P1,P2,v,F,K,llft]=llf_kfilt(y,x1,P1,H,T,V,W,z);
         KKK(:,:,im)=K; %store K,F,v as these are required by smoothing procedure
         FFF(:,:,im)=F\eye(size(F,1));
         vvv(:,:,im)=v;
         WWW(:,:,im)=W;
         PPP(:,:,im)=P2;       
     else %for non-quarterly obs, i.e. observed variable is missing
         x2=T*x1+z; %prediction of x
         P2=T*P1*T'+W; %covariance matrix of x 
         x1 = x2; 
         P1 = P2;
         WWW(:,:,im)=W;
         PPP(:,:,im)=P2;       
     end     
     
 end
 
%% DURBIN AND KOOPMAN ALGORITHM
%backward recursion for rt
T=size(Y_m,1);
rrr=zeros(4,T); %starting value set to 0

for t=T-1:-1:1
    if isnan(Y_ext(t+1,i)) %if missing quart obs for t+1
        rrr(:,t)=fff(:,:,t+1)'*rrr(:,t+1);
    else
        rrr(:,t)=(fff(:,:,t+1)-KKK(:,:,t+1)*HHH(:,t+1)')'*rrr(:,t+1)+HHH(:,t+1)*FFF(:,:,t+1)*vvv(:,t+1); %eq (11) for r(t-1) on lside 41 of Harey
    end
end
%iterate one more time to get r0
if isnan(Y_ext(1,i)) %if missing quart obs
    rrr0=fff(:,:,1)'*rrr(:,1); %if 1st obs (t+1) missing
else
    rrr0=(fff(:,:,1)-KKK(:,:,1)*HHH(:,t+1)')'*rrr(:,1)+HHH(:,t+1)*FFF(:,:,1)*vvv(:,1); %eq (11) for r(t-1) on lside 41 of Harey
end

%forward recursion for states (aaa_t) 
aaa=NaN(4,T);
a0=zeros(4,1)+zzz(:,1); %starting value of a (1 is valid if we use detrended data here)
P0=PPP(:,:,1);
r0=zeros(4,1);
aaa(:,1)=a0+P0*rrr0; %starting value for aaa
for t=2:T
        aaa(:,t)=zzz(:,t)+fff(:,:,t)*aaa(:,t-1)+sqrt(WWW(:,:,t))*sqrt(WWW(:,:,t))'*rrr(:,t-1); %%eq 13 on slide 42 of Harvey, Q=I so omitted
end        
        
%%
%multiply qt by smt (monthly trend) to get actual non-detrended monthly
%values
tmp2=aaa(1,:)'.*Y_m(:,i);

