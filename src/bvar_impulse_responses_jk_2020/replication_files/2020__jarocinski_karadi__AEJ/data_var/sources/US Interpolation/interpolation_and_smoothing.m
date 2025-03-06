function [tmp]=interpolation_and_smoothing(b_st)


   global Yreg; 
   global Xreg2;
   global Y_m;
   global Y_q;
   global i;
   global element;




tmp = max(abs(Yreg(:,i)));
 vague = 1000*tmp*tmp; 
 
 %size(Xreg2,2) stores e
 %for j=1:size(Xreg2,2) from 1 to 3
     
 small = 1.0e-08;
 large = 1.0e+02;


 h = zeros(4,1);
 r = small;
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
     %b(q).second=b_st(1:size(Xreg2(q).second,2)+2);
     %index=size(Xreg2(q).(element),2)
     b(q).(element)=b_st(index_size+1 : index_size+(size(Xreg2(q).(element),2)+2)) ;
     index_size=index_size+size(Xreg2(q).(element),2)+2;
 end
     
 
 Zreg=NaN(size(Xreg2(1).(element),1),numind);
 
 for q=1:numind
     Zreg(:,q)=Xreg2(q).(element)*[b(q).(element)(1:size(Xreg2(q).(element),2))]; %zreg1 = xreg1*bx1; each colum of Zreg contains fitted value from one set of regressors
 end
     
%  zreg1 = xreg1*bx1;
%  zreg2 = xreg2*bx2;


%%

%@ Carry Out Interpolation Using Estimated Parameters @

   
 %@ Initial Values @
 %set starting values for kalman vectors (intitially set everything to 0)
 nstate = size(f,1);
 %nm = size(xreg1,1);
 nm = size(Zreg,1); %note that xreg1 contains also NA values, hence here it doe snot matter whether we use xreg1 or xreg2 to 
 %determine the matrix size - both have equal number of rows - full 693
 p1 = vague*eye(nstate); %nstate = rows(f)=4;
 x1 = zeros(nstate,1);
 x1t=zeros(nm+1,nstate); %nstate = rows(f); so it is 666+1x4 matrix
 p1t=zeros(nm+1,nstate*nstate);
 x2t=zeros(nm+1,nstate);
 p2t=zeros(nm+1,nstate*nstate);
 x1t(1,:)=x1';  %fill the first row of x1t and p1t with starting values values
 p1t(1,:)=(p1(:))';
 
 
 
 %%
 
 
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
         x2=f*x1+z; %prediction of x
         p2=f*p1*f'+qju; %covariance matrix of x 
         x1 = x2; 
         p1 = p2;
     %NOTE: log likelihood needs to be consturcted for every period
     %separately since matrices differ across time periods (Ft and vt)
     
     %NOTE: joint density is available only for periods where y is observed
     %(since we need the updating step): hence only for months 3,6,9 and 12
     end
          
     
      x1t(im+1,:)=x1'; %store updates of x1 (predictions if no update available)
      p1t(im+1,:)=(p1(:))'; %stores updates of p1 (predictions if no update available)
      x2t(im+1,:)=x2'; %stores predictions (of unobserved component) only
      p2t(im+1,:)=(p2(:))'; %stores predictions (of covariance matrix of x) only
      
      im=im+1; 
 
 end
 

 
 %%
 
%store f matrices for both sub-samples (as based on monthly indicators)
f_base=repmat(f,numind,1);

sum=0;
for d=1:numind %create common vector for for f - i.e. put estimates of f into one matrix and the ndraw from it when&what needed
    f_base(d*4,4)= b_st(sum+size(Xreg2(d).(element),2)+1,1);
    sum=sum+size(Xreg2(d).(element),2)+2;
end

% f2 = f; 
% f2(4,4) = rho2;
 
%@ -- Kalman Smoother -- @
x3t=zeros(nm+1,nstate); %x3t and p3t use dto store smoothed estimates (of unobserved component and its covariance)
p3t=zeros(nm+1,nstate*nstate);
%fill in the last rows with most recent (last, i.e. for T) estimates of x
%and P (to start backward recursion)
x3t(nm+1,:)=x1';
p3t(nm+1,:)=(p1(:))';
x3=x1; p3=p1; %intitial values of x3 and p3 are starting vectors (i.e. those obtained for period T by KF)

j=nm; %nm is number of observations in xreg1 vector

 
%%
while j >= 1
     x2=x2t(j+1,:)'; %recover last est of x2 and transpose it to get a column vector (to ge tupdates of x2)
     p2=reshape(p2t(j+1,:),nstate,nstate); %recover est of P1, reshape needed to ge tit as a square matrix
     x1=x1t(j,:)'; %the same for predcitions 
     p1=reshape(p1t(j,:),nstate,nstate);
%          if max(isnan(zreg2(i)))>0
%          %@ Use zreg2 == invtchange, unless it is missing @
%              ft=f1;
%          else
%              ft=f2;
%          end
%          %x1, x2, p1 and p2 are estimates for last period
%          %ft is a relevant form of f matrix (depending on which sub-sample
%          %we use
            Fill=0;
            q=size(Zreg,2); %determines which index of exog regressors is considered - for the next 'while' loop
            while Fill<1
                 if max(isnan(Zreg(j,q)))==0 %if Zreg(im,q) is non-NaN - fill in the position and change Fill=1 - terminate loop/(search)
                     Fill=1;
                     ft=f_base((q-1)*4+1:4*q,:);
                 else %if Zreg(im,q) is NaN - leave Fill unchanged at 0 to ocntinue search and 
                      %reduce q by 1 (look in different column of Z i.e. try different
                      %indep regressor)
                     q=q-1;
                 end
            end
     [x3,p3]=ksmooth(x1,x2,x3,p1,p2,p3,ft);
     x3t(j,:)=x3'; %the obtained smoothed estimates of x  and p (x3 and p3) are added to vectors x3t and p3t and become input in another iteration
     p3t(j,:)=(p3(:))';
     
     if j==5
        t=0
     end
     
     j=j-1;
     
end

%%


tmp = x3t(2:size(x3t,1),1).*Y_m(:,i); %extract qt estimates from x3t (skip the first row since it contains starting values which we set and are not data based
%multiply qt by smt (monthly trend) to get actual non-detrended monthly
%values

%xlswrite('sixth indicator', tmp);
