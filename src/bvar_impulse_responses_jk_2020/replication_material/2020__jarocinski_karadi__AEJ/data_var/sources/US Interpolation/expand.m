function [Y_ext]=expand_QtoM(y, months)
%changes quarterly vector y to a monthly vector (with number of months
%specified by 'months' variable), i.e. for quarters 3,6,9,12 it assigns
%quarterly values, all other months have values 'N/A'
if size(y,1)<months %if quarterly dat aentered expand their vector to monthly size, else skip this step (as it already has the required size)
    %%
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
