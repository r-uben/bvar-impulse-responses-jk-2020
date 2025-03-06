function [Yreg, Xreg]=detrend(Y, X, Y_m, Y_q, X_m)
%function detrends the monthly and quarterly observations by trends
%specified in matrices Y_m, Y_q, X_m. Detrending is performed by dividing
%value of origina lseries by its trend
%Output: Yreg are detrended values of Y at quarterly frequency, X are detrended
%values of indep regressors at monthly frequency

%create matrices of relevant trends
Y_trend=Y_q;
X_trend=NaN(size(X,1), size(X,2));
%first indicator:
X_trend(:,2)=Y_m(:,2); X_trend(:,3)=Y_m(:,2); %second indicator: Inv in Non-residential Structures
X_trend(:,4:8)=X_m(:,4:8); %third indicator: Inv, Equipment and Software
X_trend(:,9:10)=repmat(Y_m(:,4),1,2); %fourth indicator: Residential Structures
X_trend(:,11:12)=repmat(Y_m(:,5),1,2); %fifth indicator: Change in inventories
X_trend(:,13)=Y_m(:,6); %sixth indicator (I): Exports
X_trend(:,14:16)=X_m(:,14:16); %sixth indicator (II): Exports
X_trend(:,17:18)=repmat(Y_m(:,7),1,2); %seventh indicator (I): Imports
X_trend(:,19:20)=X_m(:,19:20); %seventh indicator (II): Imports
X_trend(:,21:25)=X_m(:,21:25); %eight indicator: Government
X_trend(:,26)=X_m(:,26); %PCEPI

%Y_trend and X_trend are now matrices of trends for Y and X respectively

%detrend:
Y_trend=[zeros(size(Y_trend,1),2), Y_trend];
Yreg=Y./Y_trend;
Xreg=X./X_trend;

end