% This script prepares the environment for the estimation of shocks from 
% principal components of GSS variables
% possibly with some GSS variables included as is.

% select the variables and the number of principal components
% convention:
% X - variables from which we extract PC
% Z - variables we include as-is
% Y = [PC Z] - variables on which we estimate the model
% XZ = [X Z] - variables for which we track the effect

% Prior to running this script
% define the variable names and the number of principal components:
%
% xnames - names of X variables (string array)
% znames - names of Z variables (string array)
% Npc - number of principal components (integer)

% load GSS raw data into tab
run("../data/GSSrawdata.m")

% load the approximate maturities of the interest rates in GSS data
tab_maturities = readtable("../data/GSS_maturities.csv","ReadRowNames",true);

% prepare Z
Z = tab{:,znames};
zmaturities = tab_maturities{znames,:};

% prepare X
X = tab{:,xnames};
xmaturities = tab_maturities{xnames,:};

% pca
Xstd = std(X);
Xnorm = X./Xstd;
[coeff,score] = pca(Xnorm, 'Centered', false, 'NumComponents',Npc);
%Xrec = score*coeff'.*Xstd; % reconstructed X
%Xcentered = score*coeff';

% prepare Y
Y = [score Z];
ynames = [compose("PC%d",1:Npc) znames];
ymaturities = [(1:Npc)'; zmaturities]; % arbitrary, just so plots don't crash

% prepare info on V
xznames = [xnames, znames];
xzmaturities = [xmaturities; zmaturities];

% mapping from Y to V
Y2XZ = blkdiag(coeff'.*Xstd, eye(length(znames)));

clearvars -except Y ynames ymaturities xznames xzmaturities Y2XZ tab
