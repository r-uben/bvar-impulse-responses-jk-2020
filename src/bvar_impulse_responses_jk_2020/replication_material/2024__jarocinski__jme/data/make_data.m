% This script creates the file 'data.csv' from the source data 'GSSrawdata.xlsx'. 
% The file 'GSSrawdata.xlsx' is the dataset of Gurkaynak, Sack and Swanson (2005)
% updated until June 2019 by Gurkaynak, Karasoy-Can and Lee (2021) 
% and available in their replication files at 
% \url{http://www.bilkent.edu.tr/~refet/GKL_replication.zip}.
% Please cite this original data source.
% To run this script first download GSSrawdata.xlsx
% and save it in this folder.

% Load the updated Gurkaynak, Sack and Swanson (2005) dataset
opts = detectImportOptions('GSSrawdata.xlsx');
opts = setvartype(opts,opts.VariableNames(2:end),'double');
tab = readtable('GSSrawdata.xlsx',opts);
tab.Properties.VariableNames{1} = 'date';

% rescale to basis points
tab{:,7:end} = 100*tab{:,7:end};

% select the sample
isample = tab.date > datetime('01-Jul-1991');
isample(tab.date == datetime('17-Sep-2001')) = 0; % before mkts opened, dropped in Swanson (2021)
isample(tab.date == datetime('25-Nov-2008')) = 0; % before mkts opened, not FOMC announcement, dropped in Swanson (2021)
isample(tab.date == datetime('01-Dec-2008')) = 0; % not FOMC announcement, dropped in Swanson (2021)

tab = tab(isample,:);
fprintf('Data from %s to %s, T=%d\n', tab{1,'date'},tab{end,'date'},size(tab,1))

% select the colums and write
tab = tab(:,["date" "MP1","ONRUN2","ONRUN10","SP500"]);
tab.date.Format = 'uuuu-MM-dd';
writetable(tab, 'data.csv')