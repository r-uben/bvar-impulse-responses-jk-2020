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

global Y_m; global Y_q; global X_m;
b_start=NaN(20,size(Y_ext,2));
%First indicator: PCE
X(1).first=[];

%Second indicator: Inv in non-res structures
    % xreg1 = ones(dnobs,1)~ym_res1;
    % xreg2 = ones(dnobs,1)~ym_res2;
    X(1).second=[ones(size(Xreg,1),1), Xreg(:,2)];
    X(2).second=[ones(size(Xreg,1),1), Xreg(:,3)];

    %Starting values
    bx_start1=[0 1]';       %bx_start1 = 0|1; 
    rho_start1=0.2;         %rho_start1 = 0.2;
    tmp2 = (Y_ext(:,2)-X_tr(:,2)); % tmp = packr(nr_struc_m-pr_res1);
    tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
    seps_start1 = std(tmp)/mean(Y_m(:,2));
    
    bx_start2=[0 1]';
    rho_start2=0.2;
    tmp2 = (Y_ext(:,2)-X_tr(:,3));    % tmp = packr(nr_struc_m-pr_res2);
    tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
    seps_start2 = std(tmp)/mean(Y_m(:,2));
    
    %b_start=[];
    %b_start=[b_start, [bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2]];
    b_start1=[bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2];
    b_start(1:size(b_start1,1),2)=b_start1;
    
%     for i=1:2
%     bx_start=ones(size(X(i).second,2),1);
%     bx_start(1)=0;
%     rho_start=0.2;
%     bx_start=[bx_start; rho_start]
%     end
% 
%     tmp2 = (NR_STRUC_Q-pr_res2); 
%     tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
%     seps_start2 = std(tmp)/mean(smt); 
    
    

%Third indicator: Inv, Equipment and Software
    % xreg1 = ones(dnobs,1)~ym_atcgvs_1;
    % xreg2 = ones(dnobs,1)~ym_andevs_1~ym_aitivs_1;
    % xreg3 = ones(dnobs,1)~ym_andevs_2~ym_aitivs_2;
    X(1).third=[ones(size(Xreg,1),1), Xreg(:,4)];
    X(2).third=[ones(size(Xreg,1),1), Xreg(:,5:6)];
    X(3).third=[ones(size(Xreg,1),1), Xreg(:,7:8)];

    %Starting values
    bx_start1=[0 1]';       %@ bx are the coefficients on x @
    bx_start2=[0 1 1]';
    bx_start3=[0 1 1]';
    rho_start1=0.2;         %@ AR coefficient for errors @     
    rho_start2=0.2;
    rho_start3=0.2;
    seps_start1 = 0.3;
    seps_start2 = 0.3;
    seps_start3 = 0.3;
    
    b_start1=[bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2; bx_start3; rho_start3; seps_start3];
    b_start(1:size(b_start1,1),3)=b_start1;
    
    %b_start=[b_start, [bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2; bx_start3; rho_start3; seps_start3]];
    
  

%Fourth indicator: Residential Structures
    % xreg1 = ones(dnobs,1)~ym_confr;
    % xreg2 = ones(dnobs,1)~ym_res;
    X(1).four=[ones(size(Xreg,1),1), Xreg(:,9)];
    X(2).four=[ones(size(Xreg,1),1), Xreg(:,10)];    
    
    %Starting values
    bx_start1 = [0 1]';                             %@ bx are the coefficients on x @
    bx_start2 = [0 1]';
    rho_start1 = 0.2;                               %@ AR coefficient for errors @
    rho_start2 = 0.2;
    
    tmp1 = (Y_ext(:,4)-X_tr(:,9));    %tmp = packr(res_m-confr); 
    tmp = tmp1(~any(isnan(tmp1),2),:);  %delete rows with missing values
    seps_start1 = std(tmp)/mean(Y_m(:,4));      
    tmp1 = (Y_ext(:,4)-X_tr(:,10));    %tmp = packr(res_m-res);
    tmp = tmp1(~any(isnan(tmp1),2),:); %delete rows with missing values
    seps_start2 = std(tmp)/mean(Y_m(:,4)); 

    b_start1 = [bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2];
    b_start(1:size(b_start1,1),4)=b_start1;
    
    

%Fifth indicator: Change in Inventories
    % xreg1 = ones(dnobs,1)~ym_ivmt_ch;
    % xreg2 = ones(dnobs,1)~ym_invtchange;
    X(1).five=[ones(size(Xreg,1),1), Xreg(:,11)];
    X(2).five=[ones(size(Xreg,1),1), Xreg(:,12)];
    
    %Starting values
    bx_start1 = [0 1]';                      %@ bx are the coefficients on x @
    bx_start2 = [0 1]';
    rho_start1 = 0.2;                        %@ AR coefficient for errors @
    rho_start2 = 0.2;
    seps_start1 = 0.2;                       %@ Standard deviation of epsilon @
    seps_start2 = 0.2;                       %@ Standard deviation of epsilon @
    
    b_start1 = [bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2];
    b_start(1:size(b_start1,1),5)=b_start1;
    
    

%Sixth indicator: Exports
    % xreg1 = ones(dnobs,1)~ym_exports_bci;
    % xreg2 = ones(dnobs,1)~ym_exports_bci~ym_exports_mac~ym_exports_ag;
    % xreg3 = ones(dnobs,1)~ym_exports;
    X(1).six=[ones(size(Xreg,1),1), Xreg(:,16)];
    X(2).six=[ones(size(Xreg,1),1), Xreg(:,14:16)];
    X(3).six=[ones(size(Xreg,1),1), Xreg(:,13)];
    
    %@ Initial Values of Parameters @
    bx_start1 = [0 1]';                                %@ bx are the coefficients on x @
    bx_start2 = [0 1 1 1]';                            %@ bx are the coefficients on x @
    bx_start3 = [0 1]';
    rho_start1 = 0.2;                               %@ AR coefficient for errors @
    rho_start2 = 0.2;
    rho_start3 = 0.2;
    seps_start1 = 0.3;
    seps_start2 = 0.3;
    seps_start3 = 0.3;
    
    b_start1=[bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2; bx_start3; rho_start3; seps_start3];
    b_start(1:size(b_start1,1),6)=b_start1;
    
%Seventh indicator: Imports
    % xreg1 = ones(dnobs,1)~ym_imports_bci;
    % xreg2 = ones(dnobs,1)~ym_imports_bci~ym_imports_oil~ym_imports_auto;
    % xreg3 = ones(dnobs,1)~ym_imports;
    X(1).seven=[ones(size(Xreg,1),1), Xreg(:,20)];
    X(2).seven=[ones(size(Xreg,1),1), Xreg(:,18:20)];
    X(3).seven=[ones(size(Xreg,1),1), Xreg(:,17)];
    
    %@ Initial Values of Parameters @
    bx_start1 = [0 1]';                                %@ bx are the coefficients on x @
    bx_start2 = [0 1 1 1]';                            %@ bx are the coefficients on x @
    bx_start3 = [0 1]';
    rho_start1 = 0.2;                               %@ AR coefficient for errors @
    rho_start2 = 0.2;
    rho_start3 = 0.2;
    
    tmp1 = (Y_ext(:,7)-X_tr(:,20)); % tmp = packr(imports_m-imports_bci);
    tmp = tmp2(~any(isnan(tmp2),2),:); %delete rows with missing values
    seps_start1 = std(tmp)/mean(Y_m(:,7));
    
    seps_start2 = seps_start1;
    
    tmp1 = (Y_ext(:,7)-X_tr(:,17)); % tmp = packr(imports_m-imports);
    tmp = tmp1(~any(isnan(tmp1),2),:); %delete rows with missing values
    seps_start3 = std(tmp)/mean(Y_m(:,7));
    
    b_start1=[bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2; bx_start3; rho_start3; seps_start3];
    b_start(1:size(b_start1,1),7)=b_start1;

%Eight indicator: Government
    % xreg1 = ones(dnobs,1)~ym_wages~ym_conq;
    % xreg2 = ones(dnobs,1)~ym_wages~ym_conq~ym_man_ship_def_1;
    % xreg3 = ones(dnobs,1)~ym_wages~ym_con_gov~ym_man_ship_def_2;
    X(1).eight=[ones(size(Xreg,1),1), Xreg(:,21), Xreg(:,22)];
    X(2).eight=[ones(size(Xreg,1),1), Xreg(:,21), Xreg(:,22), Xreg(:,24)];
    X(3).eight=[ones(size(Xreg,1),1), Xreg(:,21), Xreg(:,23), Xreg(:,25)];
    
    
    %@ Initial Values of Parameters @
    bx_start1 = [0 1 1]';                         %@ bx are the coefficients on x @
    bx_start2 = [0 1 1 1]';
    bx_start3 = [0 1 1 1]';
    rho_start1 = 0.2;                            %@ AR coefficient for errors @
    rho_start2 = 0.2;
    rho_start3 = 0.2;
    seps_start1 = 0.3;
    seps_start2 = 0.3;
    seps_start3 = 0.3;

    b_start1=[bx_start1; rho_start1; seps_start1; bx_start2; rho_start2; seps_start2; bx_start3; rho_start3; seps_start3];
    b_start(1:size(b_start1,1),8)=b_start1;
    
%Ninth indicator: GDP deflator
    X(1).nine=[ones(size(Xreg,1),1), Xreg(:,26)];
    
    %@ Initial Values of Parameters @
    bx_start = [0 1]';                             %@ bx are the coefficients on x @
    rho_start = 0.2;                               %@ AR coefficient for errors @
    seps_start = 0.2;                              %@ Standard deviation of epsilon @
    b_start1 = [bx_start; rho_start; seps_start];
    b_start(1:size(b_start1,1),9)=b_start1;

    

end