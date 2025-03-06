
%%
%to count number of sheets:
% [status,sheets] = xlsfinfo(filename);
% numOfSheets = numel(sheets);

%then declasre a variable NUM=number of sheet

%for tips on looping through sheets look here:
%https://de.mathworks.com/matlabcentral/answers/57298-how-to-for-loops-for-excel-tab

%we could use some functio nto create sheet name of the sheet we want to
%open

%when a relevant sheet is opened, load it as a matrix: first series will be
%quarterly, the rest will be monthly indicators
%hence load all as vectors and work with them happily ever after

%ALTERNATIVE way: - do it This way!
%one sheet for quarterly indicators, anothe rone for monthly (this allows
%data to be updated easily)
%then create a vector which specifies how many monthly indicators
%correspond to each quarterly series and allocate them based on that (so monthly indicators need to be ordered correctly (w.r.t. quarterly indicators))

%so at first we load all the data - quarterly data into Y..., monthly into
%X..., and then based on the vector of association we draw relevant
%elements from X... for a given Y...

%There are specific data trasformations and again for easier updating, they
%should be donein MatLab
%However, this messes up the generality of script

%An option then is to create a function which, upon loading the data into X
%and Y, does all the relevant transformations, returning Y and X in the
%correct form

%once data are loaded in the transformed form ,script can be run. However,
%the number of X vairables changes the form of vectors and coefficients
%that need to be estimated by MLE - can this be easily accommodated? 

%So start with deriving the (specific) form of equations for more X
%variables (e.g. 3) and see what changes it requires throughout the code
%(as far as matrix manipulation and estimation are concerned)

%we can also have vector indicating how seps is selected, and the nloading
%the value sfro mthere (or it can conatin 1/0 value sjus tindicating the
%procedure

%not ethat no matter how many columns xreg has, prior/b-start value for the
%first is 0 and for remaining it is 1 - henc ejust identify number of
%columns and set b_start accrodingly


%% PCE
%Check whether the first monthly indicator (PCE) averages to quarterly
% PCE=X(:,1);
% PCE=reshape(PCE,[3,size(PCE,1)/3]);
% PCE=mean(PCE)';
% 
% PCE_q=Y(:,3);
% 
% Diff=PCE_q-PCE;
% S=sum(abs(Diff));
% 
% T=[PCE PCE_q abs(Diff)];
%Since it does we don't need to deal with it at all (hence all loops begin
%with the second indicator)
%%
%nknots=5;

Results=[]; %matrix storing results, i.e. smoothed and interpolated series

VoA=[1 4 5 0 2 4 4 7 1]; %specifies how many monthly indices (from Mdata) correpsond to each quarterly series (maping from Mdata to X thorugh transform)
VoA_x=[1 2 5 2 2 4 4 5 1]; %specifies correspodence of monthly indices (from X) to quarterly series (maping from X to X_m)

nknots_y=[0 5 5 5 5 4 4 5 5]; %last entries are for GDP deflator
nknots_x=[1 1 4 1 1 4 4 4 5];

[Qdata,Qtxt] = xlsread('DISTRIBUTE_GDP_GDI_INPUT2', 'Quarterly');
[Mdata,Mtxt] = xlsread('DISTRIBUTE_GDP_GDI_INPUT2', 'Monthly');

[Y, X, names]=transform(Qdata, Mdata, VoA, Qtxt, Mtxt); %function that transforms data (e.g. creat elogs, divide by 1000, first-difference...
months=size(X,1);

%compute detrended dat ausing Cubic splines
%we do separately for y and x since they have different column numbers
%(THIS MAY BE SIMPLIFIED!)
global Y_m; global Y_q; global X_m;
Y_m=[]; Y_q=[];
% for j=3:size(Y,2); %first two columns are Years and Months, so skip them
%     [Y_m(:,j-2),Y_q(:,j-2), Y_ext(:,j-2)]=CsS(Y(:,j),nknots, months); %create Y_m and Y_q - matrix of dentrended Y on monthly and quarterly frequency respectively
%     %NOTE: THIS MAY REQUIRE gaps in Y for monthly obs (e.g. months 1 2 4 5)!!!
% end

for j=3:size(Y,2); %first two columns are Years and Months, so skip them
    [Y_q(:,j-2),Y_m(:,j-2)]=cspline_qm2(Y(:,j),nknots_y(j-2)); %create Y_m and Y_q - matrix of dentrended Y on monthly and quarterly frequency respectively
    Y_ext(:,j-2)=expand(Y(:,j), months);
    %NOTE: THIS MAY REQUIRE gaps in Y for monthly obs (e.g. months 1 2 4 5)!!!
end
%fifth indicator is scaled by spline base don aboslute values of y!
W=abs(Y(:,7));
[Y_q(:,5),Y_m(:,5)]=cspline_qm2(W,nknots_y(5)); 


% nknots=5;
% for j=1:size(X,2) 
%     [X_m(:,j),X_q(:,j)]=CsS(X(:,j),nknots, months);
% end

% nknots=5;
% for j=1:size(X,2) 
%     [X_m(:,j),X_q(:,j)]=CsS(X(:,j),nknots, months);
% end

index=0; %refers to the numbe rof quarterly series t owhic X(:,j) belongs
for j=1:size(X,2)
    if j-sum(VoA_x(1:index))>0
        index=index+1;
        if VoA_x(index)==0
            index=index+1;
        end
    end
%     j
     index
     nknots_x(index)
    [X_m(:,j)]=cspline(X(:,j),nknots_x(index));
end
%sum(VoA(1:Index))
% [X_m(:,5)]=cspline(X(:,5),nknots);

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
%quarterly series

global i;
for i=2:size(b_start,2) %loop through Quarterly variables
    %i=9;
    i 
    
    global element; %element stores field name for a relevant quart/monthly indices (for reference whe nusing Xreg2)
    element=namevec{i};
    %element
    %element=namevec{8};
    
    %tmp = ones(size(Y_q,2),1);
    b_st=b_start(:,i); %load correct column vector of startng values correpsonding to the given indicator
    b_st = b_st(~any(isnan(b_st),2),:); %delete rows with missing values


    [bmax,f]=solvopt(b_st, 'lkfilt_inv_nonres5');
    
    %Five and six are problematic:
    if i==6 || i==5 || i==3 %repeat for problematic cases, taking last bmax as a starting value
        [bmax,f]=solvopt(bmax,'lkfilt_inv_nonres5');
    end
    
    %Interpolation and smoothing:
    tmp=interpolation_and_smoothing(bmax);
    
    %Save results (for components)
    Results=[Results tmp];
    
    
    
    
end

Results=[X(:,1), Results]; %add first indicator (PCE), since it was not interpolated


%% Save GDP results
%Note that here which variable is in which column matters! (sinc eimports
%needs to be substracted)

%original quarterly GDP
GDP=sum(Y(:,3:end-1),2)-2*Y(:,9); %substract imports & since they are included in the sum they must be substracted twice
GDP_m=repmat(GDP',3,1);
GDP_m=reshape(GDP_m,[size(GDP_m,2)*3,1]);

%inteprolated (monthly GDP)
%(Results and Y contain price deflator in the last column - do not include
%it in the sum)
GDP_nominal=sum(Results(:,1:end-1),2)-2*Results(:,7); %substract imports & since they are included in the sum they must be substracted twice

%Test whether interpolated quartely observations average to quarterly
%series
GDP_test=reshape(GDP_nominal,[3,size(GDP_nominal,1)/3]);
GDP_test=mean(GDP_test);
GDP_test=repmat(GDP_test,3,1);
GDP_test=reshape(GDP_test,[size(GDP_test,2)*3,1]);

GDP_real=GDP_nominal./Results(:,end)*100;

Names_Results={'GDP_m', 'GDP_test', 'GDP_nominal', 'GDP_real'}; %name Result columns
R=[GDP_m GDP_test GDP_nominal GDP_real];


 

%% TO INCORPORATE: INTERPOLATION AND SMOOTHING




%produce excel file with results
xlswrite('Interpolated_US', Results, 'Monthly Indices');
xlswrite('Interpolated_US', R, 'GDP');



