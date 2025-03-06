function [X,b_start]=regressors(Xreg, Y_ext, X_tr)
%'regressors' constructs array struct X in which all exogenous regressors
%are stored as matrices (including constant). 
%Structure of X is: X(#of regressor).#ofindicator

%NOTE: PUT regressors in increasing order: e.g. put [X(:,1),X(:,2)] before
%[X(:,1),X(:,2),X(:,3)] since the later regressors get precedence later on
%- so they should be the more informative ones
%Std_stvalue=[1 1 1 1 1 1]; %governs the function used to generate starting values for standard deviation of error

%bx_.. are starting means for monthly variables
%rho_start is an AR(1) coefficients for error in transition equation (for
%monthly indices)
%seps_start is a standard deviation of errors in AR(1) model for errors of
%trnasition equation


Std_stvalue=[1 1 1 1 1 1];

global Y_m; global Y_q; global X_m;
b_start=NaN(20,size(Y_ext,2));

%% First indicator: Private final consumption
Ind=1;

X(1).first=[ones(size(Xreg,1),1),Xreg(:,2)];
X(2).first=[ones(size(Xreg,1),1),Xreg(:,2),Xreg(:,1)];

 bx_start1 = [0 1]';                             %@ bx are the coefficients on x, i.e. (starting) means of monthly variables
 rho_start1 = 0.2;                               %@ AR coefficient for error term @
 bx_start2 = [0 1 1]';                             %@ bx are the coefficients on x @
 rho_start2 = 0.2;                               %@ AR coefficient for errors @
    
 if Std_stvalue(Ind)==1
     seps_start1= 0.3;
     seps_start2= 0.3;
 elseif Std_stvalue(Ind)==2
     tmp2 = (Y_ext(:,1)-X_tr(:,2));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start1 = std(tmp)/mean(Y_m(:,1));

     tmp2 = (Y_ext(:,1)-X_tr(:,2)-X_tr(:,1));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start2 = std(tmp)/mean(Y_m(:,1));
 elseif Std_stvalue(Ind)==3
     tmp2 = (Y_ext(:,1)-X_tr(:,2));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start1 = std(tmp)/mean(Y_m(:,1));
     
     seps_start2=seps_start1; %st.dev. of error in AR(1) model of errors of trnasition equation for monthly index
 end
     
 
    %add to b_start matrix of starting coefficients
    b_start1=[bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2];
    b_start(1:size(b_start1,1),1)=b_start1;

%% Second indicator: Government final consumption
Ind=2;

    X(1).second=[ones(size(Xreg,1),1)];

    %Starting values
    bx_start1=[1]';       
    rho_start1=0.2;         

    
 if Std_stvalue(Ind)==1
     seps_start1= 0.3;
%      seps_start2= 0.3;
%      seps_start3=0.3;
 elseif Std_stvalue(Ind)==2
     tmp2 = (Y_ext(:,2)-X_tr(:,3));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start1 = std(tmp)/mean(Y_m(:,1));

     tmp2 = (Y_ext(:,2)-X_tr(:,3)-X_tr(:,5));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start2 = std(tmp)/mean(Y_m(:,2));
     
     tmp2 = (Y_ext(:,2)-X_tr(:,3)-Xreg(:,4)-X_tr(:,5));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start3 = std(tmp)/mean(Y_m(:,2));
 elseif Std_stvalue(Ind)==3
     tmp2 = (Y_ext(:,2)-X_tr(:,3));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start1 = std(tmp)/mean(Y_m(:,2));
     
     seps_start2=seps_start1;
     seps_start3=seps_start1;
 end
    
    %add to b_start matrix of starting coefficients
    b_start1=[bx_start1; rho_start1; seps_start1];
    b_start(1:size(b_start1,1),2)=b_start1;

    
  
%% Third indicator: Gross fixed capital formation
Ind=3;
    
    X(1).third=[ones(size(Xreg,1),1), Xreg(:,7)];

    %Starting values
    bx_start1=[0 1]';       %@ bx are the coefficients on x @
    %bx_start2=[0 1 1 1]';
    rho_start1=0.2;         %@ AR coefficient for errors @     
    %rho_start2=0.2;
  
    
 if Std_stvalue(Ind)==1
     seps_start1= 0.3;
     %seps_start2= 0.3;
 elseif Std_stvalue(Ind)==2
     tmp2 = (Y_ext(:,3)-X_tr(:,7));
     %tmp2 = (Y_ext(:,3)-X_tr(:,7)-X_tr(:,8));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start1 = std(tmp)/mean(Y_m(:,3));

 elseif Std_stvalue(Ind)==3
     tmp2 = (Y_ext(:,3)-X_tr(:,7));    
     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
     seps_start1 = std(tmp)/mean(Y_m(:,3));
     
     seps_start2=seps_start1;
 end

    %add to b_start matrix of starting coefficients
    b_start1=[bx_start1; rho_start1; seps_start1];
    b_start(1:size(b_start1,1),3)=b_start1;
      

    
%% Fourth indicator: Change in business inventories
Ind=4;

    %X(1).four=[ones(size(Xreg,1),1), Xreg(:,9), Xreg(:,10)];  
    X(1).four=[ones(size(Xreg,1),1), Xreg(:,9), Xreg(:,10)]; 
    
    %Starting values
    bx_start1 = [0 1 1]';                             %@ bx are the coefficients on x @
    rho_start1 = 0.2;                               %@ AR coefficient for errors @
    
     if Std_stvalue(Ind)==1
         seps_start1= 0.3;
     elseif Std_stvalue(Ind)==2
         %tmp2 = (Y_ext(:,4)-X_tr(:,9)-X_tr(:,10));  
         tmp2 = (Y_ext(:,4)-X_tr(:,9)-X_tr(:,10));
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start1 = std(tmp)/mean(Y_m(:,4));
     elseif Std_stvalue(Ind)==3
         tmp2 = (Y_ext(:,4)-X_tr(:,9));    
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start1 = std(tmp)/mean(Y_m(:,4));
     end

    %add to b_start matrix of starting coefficients
    b_start1 = [bx_start1; rho_start1; seps_start1];
    b_start(1:size(b_start1,1),4)=b_start1;    
  
    
%% Fifth indicator: Net exports of Goods and Services
Ind=5;

    % xreg1 = ones(dnobs,1)~ym_ivmt_ch;
    % xreg2 = ones(dnobs,1)~ym_invtchange;
    X(1).five=[ones(size(Xreg,1),1), Xreg(:,11), Xreg(:,12)];
    X(2).five=[ones(size(Xreg,1),1), Xreg(:,11), Xreg(:,12), Xreg(:,13)];
    
    %Starting values
    bx_start1 = [0 1 1]';                      %@ bx are the coefficients on x @
    bx_start2 = [0 1 1 1]';
    rho_start1 = 0.2;                        %@ AR coefficient for errors @
    rho_start2 = 0.2;
    seps_start1 = 0.3;                       %@ Standard deviation of epsilon @
    seps_start2 = 0.3;                       %@ Standard deviation of epsilon @
    
     if Std_stvalue(Ind)==1
         seps_start1= 0.3;
         seps_start2= 0.3;
     elseif Std_stvalue(Ind)==2
         tmp2 = (Y_ext(:,5)-X_tr(:,11)-X_tr(:,12));    
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start1 = std(tmp)/mean(Y_m(:,5));

         tmp2 = (Y_ext(:,5)-X_tr(:,11)-X_tr(:,12)-X_tr(:,13));    
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start2 = std(tmp)/mean(Y_m(:,5));
     elseif Std_stvalue(Ind)==3
         tmp2 = (Y_ext(:,5)-X_tr(:,11));    
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start1 = std(tmp)/mean(Y_m(:,5));

         seps_start2=seps_start1;
     end
   
    %add to b_start matrix of starting coefficients
    b_start1 = [bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2];
    b_start(1:size(b_start1,1),5)=b_start1;
    
    

%% Sixth indicator: GDP deflator
Ind=6;

    % xreg1 = ones(dnobs,1)~ym_exports_bci;
    % xreg2 = ones(dnobs,1)~ym_exports_bci~ym_exports_mac~ym_exports_ag;
    % xreg3 = ones(dnobs,1)~ym_exports;
    X(1).six=[ones(size(Xreg,1),1), X_tr(:,16)];
    %X(2).six=[ones(size(Xreg,1),1), Xreg(:,14)];
    
    %@ Initial Values of Parameters @
    bx_start1 = [0 1]';                                %@ bx are the coefficients on x @
    rho_start1 = 0.2;                               %@ AR coefficient for errors @
    
    if Std_stvalue(Ind)==1
         seps_start1= 0.3;
     elseif Std_stvalue(Ind)==2
         tmp2 = (Y_ext(:,6)-X_tr(:,16));    % tmp = packr(nr_struc_m-pr_res2);
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start1 = std(tmp)/mean(Y_m(:,6));
     elseif Std_stvalue(Ind)==3
         tmp2 = (Y_ext(:,6)-X_tr(:,15));    % tmp = packr(nr_struc_m-pr_res2);
         tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
         seps_start1 = std(tmp)/mean(Y_m(:,6));

 %        seps_start2=seps_start1;
     end

    %add to b_start matrix of starting coefficients
    b_start1=[bx_start1; rho_start1; seps_start1];
    b_start(1:size(b_start1,1),6)=b_start1;
    

end