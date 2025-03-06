function tmp=cspline(q,nk)

%proc(1) = cspline(q,nk);
% /* Compute Fitted Values from Cubic Spline 
%    
%    Q = series
%    nk = number of knots 
%    
% */
 

 tmp=q;
 tmp(isnan(tmp))=1e32;
 %tmp=tmp';
 isel=(tmp~=1e32);
 trnd =(1:size(q,1))';
 trnd=selif(trnd,isel);
 qp=q(trnd);
 %qp=q(~any(isnan(q),2),:); %delete rows with missing values
 
 %qp = q[trnd]; /* select elements of q for which */
 n=size(qp,1); %n = rows(qp);
 tr = (1:n)';
 nks = zeros(nk+1,1); %/* nk=number of knots */
 i=1; 
 while i <= nk
 	nks(i+1)=floor((i/(nk+1))*n); %indices of knot points, first row in nks has 0 in it
    i=i+1; 
 end
 
 spreg = zeros(n,nk+4);
 t0 = ones(n,1);
 t1 = tr;
 t2 = tr.^2;
 spreg(:,1) = t0;
 spreg(:,2) = t1;
 spreg(:,3) = t2;
 i = 1; 
 while i <= nk+1
	tmp = t1-nks(i); %index of current obs from knot observation
	spreg(:,3+i) = (0.5*(tmp.^3)).*(t1 > nks(i)); %columns 4,5,6,7,8,9 - distance from 0 and five knot points
    i=i+1; 
 end;
 b=(spreg'*spreg)\spreg'*qp; %b = qp./spreg;
 ct = spreg*b;
 tmp = NaN(size(q,1),1);
 tmp(trnd) = ct;
 %retp(tmp);
 end