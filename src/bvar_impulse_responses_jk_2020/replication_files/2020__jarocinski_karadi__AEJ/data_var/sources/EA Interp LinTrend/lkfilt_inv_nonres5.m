function [llfs]=lkfilt_inv_nonres5(b_st)
%consturct log-likelihood, given vector of parameters 'b_st'
%uses global variables as inputs
% The system is:
%   y_t=H*x_t+e_t
%   x_t=f*x_t-1+z+v_t
%
%   Var(e_t)=V=0
%   Var(v_t)=W

 global Yreg; global Xreg2;  global Y_m; global Y_q; global i; global element;

 b_st = b_st(~any(isnan(b_st),2),:); %delete rows with missing values
 
 llf = zeros(floor(size(Yreg(:,i),1)),1);  %yq are scaled quarterly observations (of a variable to be detrended) by quarterly trend
 tmp = max(abs(Yreg(:,i)));
 vague = 1000*tmp*tmp; 
    
 H = zeros(4,1);
 V = 0;

 T = zeros(4,4);
 T(1,4) = 1;
 T(2:3,1:2) = eye(2);
 
 %q is covariance matrix in the transition eq
 W = zeros(4,4);
 %zreg are fitted values of the regression element in the transition equation
 
 %% count number of monthly indicators from Xreg2 for a given quarterly series - store their number in numind
 numind=0;
 for j=1:size(Xreg2,2)
    if size(Xreg2(j).(element),2)>0
        numind=numind+1;
    end     
 end
 %%
 
 %Split b into components for each xreg (monthly index) separately
 index_size=0;

 b(1).(element)=NaN(1,size(Xreg2(1).(element),2)+2);
 b(2).(element)=NaN(1,size(Xreg2(1).(element),2)+2);
 
 for q=1:numind
     %q
     %b(q).second=b_st(1:size(Xreg2(q).second,2)+2);
     %index=size(Xreg2(q).(element),2)
     b(q).(element)=b_st(index_size+1 : index_size+(size(Xreg2(q).(element),2)+2)) ;
     index_size=index_size+size(Xreg2(q).(element),2)+2;
 end
     
 Zreg=NaN(size(Xreg2(1).(element),1),numind);
 for q=1:numind
     Zreg(:,q)=Xreg2(q).(element)*[b(q).(element)(1:size(Xreg2(q).(element),2))]; %zreg1 = xreg1*bx1;
 end    
  
 %@ Initial Values @
 P1 = vague*eye(size(T,1));
 x1 = zeros(size(T,1),1); %starting value for x (unobserved component) 

 %Now: fill the z matrix, i.e. the regression component in the transition
 %equation (using zreg1 and zreg2 from above)
 %zreg2 (based on xreg2 takes precedence - only if it is missing is zreg1 and
 %xreg1 used
 
 for im=1:size(Zreg,1) %i.e. while im is less or equal to the number of monthly observations     
     iq = im/3; %covnerts im into quarters: since we start in Januari, Q1 this is valid
     %iq is integer for months 3 6 9 and 12, otherwise its a fraction
     Z=NaN(4,1); %4x1 matrix of missing values 
     Fill=0; %indicates whether the position for obs im has yet been filled in by fitted values of indep regressors

     %Now fill in z matrix - external component in the transition equation
     q=size(Zreg,2); %determines which index of exog regressors is considered - for the next 'while' loop
     while Fill<1 
         if max(isnan(Zreg(im,q)))==0 %if Zreg(im,q) is non-NaN - fill in the position and change Fill=1 - terminate loop/(search)
             Z = vertcat(Zreg(im,q),zeros(3,1)); 
             T(4,4) = b(q).(element)(size(Xreg2(q).(element),2)+1);
             W(4,4) = b(q).(element)(size(Xreg2(q).(element),2)+2)^2;
             Fill=1;
         else %if Zreg(im,q) is NaN - leave Fill unchanged at 0 to ocntinue search and 
              %reduce q by 1 (look in different column of Z i.e. try different
              %indep regressor)
             q=q-1;
         end
     end 
 
     %For quarterly data points, compute log-likelihood
     if iq - floor(iq)== 0 %if iq is integer, i.e. if month=3,6,9 or 12
         y = Yreg(iq,i); %no change in y
         H(1) = Y_m(im,i)/(3*Y_q(iq,i)); %smt is a monthly trend from cubic spline
         H(2) = Y_m(im-1,i)/(3*Y_q(iq,i)); %smt and sqt re-weight to take trends into account, 1/3 come from averaging - 1 quarter has 3 months
         H(3) = Y_m(im-2,i)/(3*Y_q(iq,i));
         
         %feed data to llf function to get likelihood function for this
         %time period
         [x1,x2,P1,P2,e,F,K,llft]=llf_kfilt(y,x1,P1,H,T,V,W,Z);
         %[x1,~,p1,~,~,~,~,llft]=llf_kfilt(y,x1,p1,h,f,r,qju,z);
         llf(iq) = llft;
                  
     else
         x1=T*x1+Z; %prediction of x
         P1=T*P1*T'+W; %covariance matrix of x 
     %NOTE: log likelihood needs to be consturcted for every period
     %separately since matrices differ across time periods (Ft and vt)
     
     %NOTE: joint density is available only for periods where y is observed
     %(since we need the updating step): hence only for months 3,6,9 and 12
     end

 end

 llfs=sum(llf); %log likelihood is a sum of log-lik for all the periods
     
 
end
