function [X,b_start]=regressors(Xreg, Y_ext, X_tr)
%'regressors' constructs array struct X in which all exogenous regressors
%are stored as matrices (including constant). 
%Structure of X is: X(#of regressor).#ofindicator

%NOTE: PUT regressors in increasing order: e.g. put [X(:,1),X(:,2)] before
%[X(:,1),X(:,2),X(:,3)] since the latter regressors get precedence later on
%- so they should be the more informative ones

%bx_.. are starting means for monthly variables
%rho_start is an AR(1) coefficients for error in transition equation (for
%monthly indices)
%sigma_start is a standard deviation of errors in AR(1) model for errors of
%the transition equation (i.e. standard deviation of epsilon)


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
    
sigma_start1= 0.3;
sigma_start2= 0.3;

%add to b_start matrix of starting coefficients
b_start1=[bx_start1; rho_start1; sigma_start1; bx_start2; rho_start2; sigma_start2];
b_start(1:size(b_start1,1),1)=b_start1;

%% Second indicator: Government final consumption
Ind=2;
X(1).second=[ones(size(Xreg,1),1)];

%Starting values
bx_start1=[1]';       
rho_start1=0.2;         

sigma_start1= 0.3;
    
%add to b_start matrix of starting coefficients
b_start1=[bx_start1; rho_start1; sigma_start1];
b_start(1:size(b_start1,1),2)=b_start1;

    
  
%% Third indicator: Gross fixed capital formation
Ind=3;
    
X(1).third=[ones(size(Xreg,1),1), Xreg(:,7)];

%Starting values
bx_start1=[0 1]';       %@ bx are the coefficients on x @
rho_start1=0.2;         %@ AR coefficient for errors @     

sigma_start1= 0.3;

%add to b_start matrix of starting coefficients
b_start1=[bx_start1; rho_start1; sigma_start1];
b_start(1:size(b_start1,1),3)=b_start1;
      

    
%% Fourth indicator: Change in business inventories
Ind=4;

X(1).four=[ones(size(Xreg,1),1), Xreg(:,9), Xreg(:,10)]; 
    
%Starting values
bx_start1 = [0 1 1]';                            %@ bx are the coefficients on x @
rho_start1 = 0.2;                               %@ AR coefficient for errors @
sigma_start1= 0.3;
     
%add to b_start matrix of starting coefficients
b_start1 = [bx_start1; rho_start1; sigma_start1];
b_start(1:size(b_start1,1),4)=b_start1;    
  
    
%% Fifth indicator: Net exports of Goods and Services
Ind=5;

X(1).five=[ones(size(Xreg,1),1), Xreg(:,11), Xreg(:,12)];
X(2).five=[ones(size(Xreg,1),1), Xreg(:,11), Xreg(:,12), Xreg(:,13)];
    
%Starting values
bx_start1 = [0 1 1]';                      %@ bx are the coefficients on x @
bx_start2 = [0 1 1 1]';
rho_start1 = 0.2;                         %@ AR coefficient for errors @
rho_start2 = 0.2;
sigma_start1 = 0.3;                       %@ Standard deviation of epsilon @
sigma_start2 = 0.3;                       %@ Standard deviation of epsilon @
    
%add to b_start matrix of starting coefficients
b_start1 = [bx_start1; rho_start1; sigma_start1; bx_start2; rho_start2; sigma_start2];
b_start(1:size(b_start1,1),5)=b_start1;
    
    

%% Sixth indicator: GDP deflator
Ind=6;

X(1).six=[ones(size(Xreg,1),1), Xreg(:,15), X_tr(:,16)];
    
%Initial Values of Parameters
bx_start1 = [0 1 1]';                       %bx are the coefficients on x 
rho_start1 = 0.2;                           % AR coefficient for errors
sigma_start1= 0.3;                          % st.deviation for error in AR(1)             
    
%add to b_start matrix of starting coefficients
b_start1=[bx_start1; rho_start1; sigma_start1];
b_start(1:size(b_start1,1),6)=b_start1;
    

end