% select the variables and the number of principal components
% convention:
% X - variables from which we extract PC
% Z - variables we include as-is
% Y = [PC Z] - variables on which we estimate the model
% XZ = [X Z] - variables for which we track the effect
xnames = ["MP1","FF3","ED2","ED3","ED4","ONRUN2","ONRUN5","ONRUN10","SP500"];
Npc = 5;
znames = [];

% call the script
prepare_data

% fix the signs and ordering of the shocks
sorder = [1 -2 4 3 5];
