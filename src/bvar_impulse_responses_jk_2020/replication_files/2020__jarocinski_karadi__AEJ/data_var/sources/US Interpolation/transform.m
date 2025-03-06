function [Y, X, names]=transform(Qdata, Mdata, VoA, Qtxt, Mtxt)
%'transform' transforms given data to the specific variables needed for KF interpolations
%The function of 'transform'is specific - given new data, it should be
%revised as it generally contains commandse specific to each
%quarterly/monthly observation which cannot be incorporated into loops
Y=Qdata; %no change i ndependent variables
X=[]; %build up X - matrix of vectors of regression variables
names=[];
%% First Indicator Personal Consumption Expenditures

%PCE

Index=1-1;

 xreg1=Mdata(:,2+1+sum(VoA(1:Index)));
 X=[X, xreg1];

%str = "PCEPI"; isel = namevec .$== str;  p_pce=selif(datamat',isel)';

names=[names, 'PCE1'];

%% Second Indicator - Investment: Nonresidential Structures
Index=2-1;

xreg1 = (Mdata(:,6)-Mdata(:,4))/1000; %pr_res1 = (CONP-CONFR)/1000;
X=[X, xreg1];
xreg2 = (Mdata(:,7)-Mdata(:,5))/1000; %pr_res2 = (PRIV-RES)/1000;
X=[X, xreg2];

names=[names, ' I_NS1' ' I_NS2'];

%we drop trimming x vectors and will follow the order of precedence. This
%will be set up by order of elements in X: the first is least important,
%while the last one the most (i.e. it takes precedence over all other
%elements)

%% Third indicator Investment: Equipment and Software
Index=3-1;

 xreg1=12*(Mdata(:,2+1+sum(VoA(1:Index))))/1000;          % atcgvs_1 = 12*atcgvs_1/1000;
 xreg2=12*(Mdata(:,2+1+sum(VoA(1:Index))+1))/1000;                   % andevs_1 = 12*andevs_1/1000;
 xreg3=12*(Mdata(:,2+1+sum(VoA(1:Index))+2))/1000;                   % aitivs_1 = 12*aitivs_1/1000;
 xreg4=12*(Mdata(:,2+1+sum(VoA(1:Index))+3))/1000;                   % andevs_2 = 12*andevs_2/1000;
 xreg5=12*(Mdata(:,2+1+sum(VoA(1:Index))+4))/1000;                   % aitivs_2 = 12*aitivs_2/1000;
 xreg=[xreg1 xreg2 xreg3 xreg4 xreg5];
 
 names=[names, ' I_ES1', ' I_ES2', ' I_ES3', ' I_ES4', ' I_ES5'];
 

%OR:
%create vecro of functions and then loop thourg index of quart obs, loading
%the relevant monthly indicators according to VoA and tranforming them
%using relevant function in the function vector
 
 X=[X, xreg];
     
 %% Fourth indicator: Investment in Residential Structures
 Index=4-1;

xreg1=Mdata(:,4)/1000;          %  confr = confr/1000;        @ Billions at an annual rate @ 
xreg2=Mdata(:,5)/1000;          %  res = res/1000;            @ Billions at an annual rate @

xreg=[xreg1 xreg2];

names=[names, ' I_RS1', ' I_RS2'];

X=[X, xreg];
 
 %% Fifth indicator Investment: Change in private inventories
 Index=5-1;
 xreg1=Mdata(:,2+1+sum(VoA(1:Index)))/1000;           % invtchange = inv_ch/1000;
 xreg2=Mdata(:,2+1+sum(VoA(1:Index))+1);
 xreg2=vertcat(0.041,[xreg2(2:size(xreg2,1))]-xreg2(1:(size(xreg2,1)-1)));           % ivmt_ch = vertcat(0.041,(ivmt(2:size(ivmt,2))]-ivmt[1:(size(ivmt,1)-1)]));
 xreg2=12*xreg2;                   % ivmt_ch = 12*ivmt_ch;  %@ At annual rate @

 xreg=[xreg1 xreg2];
 
 names=[names, ' I_chPI1', ' I_chPI2'];
 
 X=[X, xreg];
            
 %% Sixth Indicator: Exports
 Index=6-1;
 
 
 xreg1=12*Mdata(:,(2+1+sum(VoA(1:Index))))/1000;            %exports = 12*exports/1000;  %@ Billions At annual rate @
 xreg2=12*Mdata(:,(2+1+sum(VoA(1:Index))+1))/1000;           %exports_mac = 12*exports_mac/1000;  %@ Billions At annual rate @
 xreg3=12*Mdata(:,(2+1+sum(VoA(1:Index))+2))/1000;           %exports_ag = 12*exports_ag/1000;  %@ Billions At annual rate @
 xreg4=12*Mdata(:,(2+1+sum(VoA(1:Index))+3))/1000;           %exports_bci = 12*exports_bci/1000;  %@ Billions At annual rate @           
            
 xreg=[xreg1 xreg2 xreg3 xreg4];
 
 names=[names, ' X1', ' X2', ' X3', ' X4'];
 
 X=[X xreg];
 
 %% Seventh indicator: Imports
 Index=7-1;
 

xreg1=12*Mdata(:,(2+1+sum(VoA(1:Index))))/1000;            % imports = 12*imports/1000;  %@ Billions At annual rate @
xreg2=12*Mdata(:,(2+1+sum(VoA(1:Index))+1))/1000;            % imports_oil = 12*imports_oil/1000;  %@ Billions At annual rate @
xreg3=12*Mdata(:,(2+1+sum(VoA(1:Index))+2))/1000;            % imports_auto = 12*imports_auto/1000;  %@ Billions At annual rate @
xreg4=12*Mdata(:,(2+1+sum(VoA(1:Index))+3))/1000;            % imports_bci = 12*imports_bci/1000;  %@ Billions At annual rate @
 
xreg=[xreg1 xreg2 xreg3 xreg4];

names=[names, ' IM1', ' IM2', ' IM3', ' IM4'];

X=[X xreg];

%% Eight Indicator: Government
Index=8-1;

xreg1=Mdata(:,(2+1+sum(VoA(1:Index))));              % Wages_G
xreg2=Mdata(:,(2+1+sum(VoA(1:Index))+1))/1000;       % conq = conq/1000;
xreg3=Mdata(:,(2+1+sum(VoA(1:Index))+2))/1000;       % con_gov = con_gov/1000;
xreg4=Mdata(:,(2+1+sum(VoA(1:Index))+3))-Mdata(:,(2+1+sum(VoA(1:Index))+4));      % man_ship_def_1 = amtmvs_1 - amxdvs_1;
xreg5=Mdata(:,(2+1+sum(VoA(1:Index))+5))-Mdata(:,(2+1+sum(VoA(1:Index))+6));      % man_ship_def_2 = amtmvs_2 - amxdvs_2;
xreg4=12*xreg4/1000;            % man_ship_def_1 = 12*man_ship_def_1/1000;
xreg5=12*xreg5/1000;            % man_ship_def_2 = 12*man_ship_def_2/1000;


xreg=[xreg1 xreg2 xreg3 xreg4 xreg5];

names=[names, ' G1', ' G2', ' G3', ' G4', ' G5'];

X=[X xreg];

%% Ninth Indicator: ???
Index=9-1;

%str = "PCEPI"; isel = namevec .$== str;  p_pce=selif(datamat',isel)';

xreg1=Mdata(:,(2+1+sum(VoA(1:Index))));            %- no ln sinc elog not used in estimation
%xreg1=log(Mdata(:,(2+1+sum(VoA(1:Index)))));            %ln_p_pce = ln(p_pce);

X=[X xreg1];

names=[names, ' nine'];

%% Quarterly Series

Y=[];

%%
Y=Qdata;

end
