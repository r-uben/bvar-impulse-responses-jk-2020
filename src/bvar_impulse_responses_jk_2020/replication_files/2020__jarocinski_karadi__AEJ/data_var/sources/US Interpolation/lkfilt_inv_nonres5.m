function [llfs]=lkfilt_inv_nonres4(b_st)
%proc lkfilt_inv_nonres(b,ytmp);
%b
   global Yreg; 
   global Xreg2;
   global Y_m;
   global Y_q;
   global i;
   global element;
 
 b_st = b_st(~any(isnan(b_st),2),:); %delete rows with missing values
 
 llf = zeros(floor(size(Yreg(:,i),1)),1);  %yq are scaled quarterly observations (of a variable to be detrended) by quarterly trend
 tmp = max(abs(Yreg(:,i)));
 vague = 1000*tmp*tmp; 
 
 h = zeros(4,1);
 r = 0;
 %create f for transition matrix,which has a form
%  0 0 0 1
%  1 0 0 0
%  0 1 0 0
%  0 0 0 rho
 f = zeros(4,4);
 f(1,4) = 1;
 f(2:3,1:2) = eye(2);
 
 %q is covariance matrix in the transition eq
 qju = zeros(4,4);
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
 p1 = vague*eye(size(f,1));
 x1 = zeros(size(f,1),1); %starting value for x (unobserved component) 
 im=1; 
 %Now: fill the z matrix, i.e. the regression component in the transition
 %equation (using zreg1 and zreg2 from above)
 %zreg2 (based on xreg2 takes precedence - only if it is missing is zreg1 and
 %xreg1 used
 
 while im<=size(Zreg,1) %i.e. while im is less or equal to the number of monthly observations
     
     iq = im/3; %covnerts im into quarters: since we start in Januari, Q1 this is valid
     %iq is integer for months 3 6 9 and 12, otherwise its a fraction
     z=NaN(4,1); %4x1 matrix of missing values 
     Fill=0; %indicates whether the position for obs im has yet been filled in by fitted values of indep regressors

     %Now fill in z matrix - external component in the transition equation
     q=size(Zreg,2); %determines which index of exog regressors is considered - for the next 'while' loop
     while Fill<1 
         if max(isnan(Zreg(im,q)))==0 %if Zreg(im,q) is non-NaN - fill in the position and change Fill=1 - terminate loop/(search)
             z = vertcat(Zreg(im,q),zeros(3,1)); 
             f(4,4) = b(q).(element)(size(Xreg2(q).(element),2)+1);
             qju(4,4) = b(q).(element)(size(Xreg2(q).(element),2)+2)^2;
             Fill=1;
         else %if Zreg(im,q) is NaN - leave Fill unchanged at 0 to ocntinue search and 
              %reduce q by 1 (look in different column of Z i.e. try different
              %indep regressor)
             q=q-1;
         end
     end
     
     %If there is a bug and z is empty, return error:
     if max(isnan(z))>0
         %"Missing data in z ... stopping";;stop;
         fprintf('Missing data in z ... stopping \n');
         quit cancel;    
     end     
 
     %For quarterly data points, compute log-likelihood
     if iq - floor(iq)== 0 %if iq is integer, i.e. if month=3,6,9 or 12
         y = Yreg(iq,i); %no change in y
         h(1) = Y_m(im,i)/(3*Y_q(iq,i)); %smt is a monthly trend from cubic spline
         h(2) = Y_m(im-1,i)/(3*Y_q(iq,i)); %smt and sqt re-weight to take trends into account, 1/3 come from averaging - 1 quarter has 3 months
         h(3) = Y_m(im-2,i)/(3*Y_q(iq,i));
         
         %feed data to llf function to get likelihood function for this
         %time period
         [x1,x2,p1,p2,e,F,llft]=llf_kfilt(y,x1,p1,h,f,r,qju,z);
         llf(iq) = llft;
                  
     else
         x1=f*x1+z; %prediction of x
         p1=f*p1*f'+qju; %covariance matrix of x 
     %NOTE: log likelihood needs to be consturcted for every period
     %separately since matrices differ across time periods (Ft and vt)
     
     %NOTE: joint density is available only for periods where y is observed
     %(since we need the updating step): hence only for months 3,6,9 and 12
     end
         im=im+1; 

 end

 llfs=sum(llf); %log likelihood is a sum of log-lik for all the periods


end
