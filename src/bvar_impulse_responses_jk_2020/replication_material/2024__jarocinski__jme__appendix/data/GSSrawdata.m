% Load the updated Gurkaynak, Sack and Swanson (2005) dataset, 
% from http://www.bilkent.edu.tr/~refet/GKL_replication.zip
opts = detectImportOptions('GSSrawdata.xlsx');
opts = setvartype(opts,opts.VariableNames(2:end),'double');
tab = readtable('GSSrawdata.xlsx',opts);
tab.Properties.VariableNames{1} = 'date';

% rescale to basis points
tab{:,7:end} = 100*tab{:,7:end};

% select the sample
isample = tab.date > datetime('01-Jul-1991');
isample(tab.date == datetime('17-Sep-2001')) = 0; % before mkts opened, dropped in Swanson (2020)
isample(tab.date == datetime('25-Nov-2008')) = 0; % before mkts opened, not FOMC announcement, dropped in Swanson (2020)
isample(tab.date == datetime('01-Dec-2008')) = 0; % not FOMC announcement, dropped in Swanson (2020)

tab = tab(isample,:);
fprintf('Data from %s to %s, T=%d\n', tab{1,'date'},tab{end,'date'},size(tab,1))

clear opts isample