function [Y_ext]=expand_QtoM(y, months)
%x=expand(y,months)
%changes quarterly vector y to a monthly vector x (with number of months
%specified by 'months' variable), i.e. for quarters 3,6,9,12 it assigns
%quarterly values, all other months have values 'N/A'
if size(y,1)<months %if quarterly dat aentered expand their vector to monthly size, else skip this step (as it already has the required size)
    %%
    %we will now try to start with quart vector - we need to creat emonthly
    %vector out of it (with months 1,2,4,6,7,8 empty)
    % quarters=repmat([3 6 9 12]',floor(size(y,1)/4)+1,1);
    % quarters=quarters(1:size(y,1),1); %trim quarters so that it is the same length as y
    % y=[y,quarters]; %first column of y now contains y, while second has quarterly positions containinhg months (3, 6, 9, 12 etc)

    %we need to import monthly size: "months"
    Z=(1:months)';
    Z=[Z,zeros(months,1)];

    j=1; %tracks quarters in y
    for i=1:size(Z,1)
        if mod(Z(i,1),3)~=0 %if month other than 3 6 9 or 12
            Z(i,2)=NaN;
        else %if end-quarter month - fill in from y
            Z(i,2)=y(j,1);
            j=j+1; %update quarterly index
        end
    end

    y=Z(:,2); %second column contains quart values (and NAs), first codes of months
    Y_ext=Z(:,2); %Y_ext gives extended vectors of quarterly data (with NaNs for months 1 2 4 5 ...)
end

end
