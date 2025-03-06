%% This code carries out an interpolation of EA GDP fro mquaterly to monthly values
% detrending of series is done using linear trend


% Load data from Excel
[Y,Qtxt] = xlsread('Input', 'Quarterly');
[X,Mtxt] = xlsread('Input', 'Monthly');

months=size(X,1);

%%
%compute detrended data using Cubic splines
%we do separately for y and x since they have different column numbers

global Y_m; global Y_q; global X_m; global Y_ext;
Y_m=[]; Y_q=[]; %Y_m and Y_q store monthly and quarterly trend respectively, extracted by cubic splines

for j=1:size(Y,2); 
    [Y_q(:,j),Y_m(:,j)]=LinTrend_qm(Y(:,j)); %create Y_m and Y_q - matrix of detrended Y on monthly and quarterly frequency respectively
    Y_ext(:,j)=expand(Y(:,j), months); %Y_ext is a monthly vector contianing only quarterly values, 
    %i.e. months 1 2 4 5 etc have missing values, while 3 6 9 12 have quarterly values
    
end

%Compute trends for X:
index=0; %refers to the number of quarterly series to which X(:,j) belongs
for j=1:size(X,2)
    [X_m(:,j)]=LinTrend(X(:,j));
end

global Yreg; global Xreg;
%now detrend the variables as appropriate, trends come from Y_m, Y_q, X_m
Yreg=Y./Y_q;
Xreg=X./X_m;

global Xreg2;
%[Xreg2,b_start]=regressors_EA_onlyHICP(Xreg, Y_ext, X); %put regressors together and get starting values (for MLE optimization)
[Xreg2,b_start]=regressors_EA2(Xreg, Y_ext, X); %put regressors together and get starting values (for MLE optimization)

%% Optimization to recover KF parameters by MLE, Inteprolation and Smoothing
namevec={'first','second','third','four','five','six','seven','eight','nine'}; %vector of names of elements of Xreg2 - needed
%in order to loop through its fields (=sets of monthly indicators for quarterly series)

Results=NaN(size(Y_ext,1),size(Y_ext,2)); %matrix storing results, i.e. smoothed and interpolated series
global i;
for i=1:size(b_start,2) %loop through Quarterly variables
    i
    global element; %element stores field name for a relevant quart/monthly indices (for reference whe nusing Xreg2)
    element=namevec{i};

    b_st=b_start(:,i); %load correct column vector of startng values correpsonding to the given indicator
    b_st = b_st(~any(isnan(b_st),2),:); %delete rows with missing values


    [bmax,f]=solvopt(b_st, 'lkfilt_inv_nonres5');
    
    %Repeat for problematic cases, using resulting bmax as starting values
    %(this approach follows suggestion by solvopt command)
    if i==5 %repeat for problematic cases, taking last bmax as a starting value
        [bmax,f]=solvopt(bmax,'lkfilt_inv_nonres5');
    end
    bmax
    
    %Interpolation and smoothing:
    [tmp2]=interpolation_and_smoothing_third2(bmax);
    
    %Save results (for components)
    Results(:,i)=[tmp2];
        
end


%% Save GDP results
%This section differs based on what quarterly indicators we use, since if
%imports are used separately we must substract them. For EA we use net
%export balance - no deductions needed

%Inteprolated (monthly GDP)
%(Results contain price deflator in the last column - do not include it in the sum)
GDP_nominal=sum(Results(:,1:end-1),2); %nominal GDP is a sum of quarterly components (excluding GDP deflator)


%Inteprolated real GDP
GDP_real=GDP_nominal./Results(:,end)*100;

%Results
R=[GDP_nominal GDP_real Results(:,end)];

%% Plot original nominal GDP and nominal GDP
plot(GDP_nominal)

%% Export results to excel 
xlswrite('Interpolated_EA_LinTrend', R, 'GDP_interpolated'); % exports interpolated nGDP, rGDP and GDP deflator
xlswrite('Interpolated_EA_LinTrend', Results, 'Interpolated_Components'); % exports imnterpolated monthly components + GDP deflator
