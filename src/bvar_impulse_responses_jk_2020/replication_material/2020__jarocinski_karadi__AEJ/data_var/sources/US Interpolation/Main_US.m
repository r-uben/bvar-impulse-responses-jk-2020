%% This code carries out an interpolation of US GDP fro mquaterly to monthly values
% detrending of series is done by cubic splines
% Approach follows the one described by Stock & Watson
% Note that for the first monthly category (PCE) monthly data average to
% quarterly - no need for interpolation

%% PCE - 
%Check whether the first monthly indicator (PCE) averages to quarterly

%%
Results=[]; %matrix storing results, i.e. smoothed and interpolated series

VoA=[1 4 5 0 2 4 4 7 1]; %specifies how many monthly indices (from Mdata) correpsond to each quarterly series (maping from Mdata to X thorugh transform)
VoA_x=[1 2 5 2 2 4 4 5 1]; %specifies correspodence of monthly indices (from X) to quarterly series (maping from X to X_m)

nknots_y=[0 5 5 5 5 4 4 5 5]; %last entries are for GDP deflator
nknots_x=[1 1 4 1 1 4 4 4 5];

%Load data
[Qdata,Qtxt] = xlsread('DISTRIBUTE_GDP_GDI_INPUT2', 'Quarterly');
[Mdata,Mtxt] = xlsread('DISTRIBUTE_GDP_GDI_INPUT2', 'Monthly');

[Y, X, names]=transform(Qdata, Mdata, VoA, Qtxt, Mtxt); %function that transforms data (e.g. creat elogs, divide by 1000, first-difference...
months=size(X,1);

%compute detrended dat ausing Cubic splines
%we do separately for y and x since they have different column numbers
global Y_m; global Y_q; global X_m;
Y_m=[]; Y_q=[];

for j=3:size(Y,2); %first two columns are Years and Months, so skip them
    [Y_q(:,j-2),Y_m(:,j-2)]=cspline_qm2(Y(:,j),nknots_y(j-2)); %create Y_m and Y_q - matrix of dentrended Y on monthly and quarterly frequency respectively
    Y_ext(:,j-2)=expand(Y(:,j), months);
end
%fifth indicator is scaled by spline based on aboslute values of y!
W=abs(Y(:,7));
[Y_q(:,5),Y_m(:,5)]=cspline_qm2(W,nknots_y(5)); 

index=0; %refers to the number of quarterly series to which X(:,j) belongs
for j=1:size(X,2)
    if j-sum(VoA_x(1:index))>0
        index=index+1;
        if VoA_x(index)==0
            index=index+1;
        end
    end
    [X_m(:,j)]=cspline(X(:,j),nknots_x(index));
end

global Yreg; global Xreg;
%now detrend the variables as appropriate
[Yreg, Xreg]=detrend(Y, X, Y_m, Y_q, X_m); 
Yreg=Yreg(:,3:end);

global Xreg2;
[Xreg2,b_start]=regressors(Xreg, Y_ext, X); %put regressors together and get starting values (for MLE optimization)

global Y_data; %drop the year and month from Y - Y_data are 223 obs of quarterly values of Y
Y_data=Y(:,3:end);

%% Optimization to recover KF parameters by MLE, Inteprolation and Smoothing

namevec={'first','second','third','four','five','six','seven','eight','nine'}; %vector of names of elements of Xreg2 - needed
%in order to loop through its fields (=sets of monthly indicators for
%quarterly series)

global i;
for i=2:size(b_start,2) %loop through Quarterly variables
    i 
    
    global element; %element stores field name for a relevant quart/monthly indices (for reference whe nusing Xreg2)
    element=namevec{i};
    
    b_st=b_start(:,i); %load correct column vector of startng values correpsonding to the given indicator
    b_st = b_st(~any(isnan(b_st),2),:); %delete rows with missing values

    [bmax,f]=solvopt(b_st, 'lkfilt_inv_nonres5');
    
    %Repeat for problematic cases
    if i==6 || i==5 || i==3 %repeat for problematic cases, taking last bmax as a starting value
        [bmax,f]=solvopt(bmax,'lkfilt_inv_nonres5');
    end
    
    %Interpolation and smoothing:
    tmp=interpolation_and_smoothing(bmax);
    
    %Save results (for monthly components)
    Results=[Results tmp];
 
end

Results=[X(:,1), Results]; %add first indicator (PCE), since it was not interpolated

%% Save GDP results
%Note that here which variable is in which column matters! (sinc eimports
%needs to be substracted)

GDP_nominal=sum(Results(:,1:end-1),2)-2*Results(:,7); %substract imports & since they are included in the sum they must be substracted twice
GDP_real=GDP_nominal./Results(:,end)*100;

R=[GDP_nominal GDP_real Results(:,end)];


%% Export results to excel file
xlswrite('Interpolated_US', Results, 'Monthly Indices');
xlswrite('Interpolated_US', R, 'GDP');



