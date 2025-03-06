function [sqt, smt]=cspline_qm2(q,nk);

 n = size(q,1); 
 tr = (1:n);
 nks = zeros(nk+1,1);
 i=1; 
 while i <= nk;
 	nks(i+1)=floor((i/(nk+1))*n);
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
 while i <= nk+1;
	tmp = t1-nks(i);
	spreg(:,3+i) = (0.5*(tmp.^3)).*(t1 > nks(i));
    i=i+1; 
 end
 %b = q/spreg;
 b=inv(spreg'*spreg)*spreg'*q;
 %print(spreg);
 sqt = spreg*b;

 trm = (1:3*n)/3;
 trm = trm + 1/3;              %@ Puts quarterly value in middle of quarter @
 spreg_m = zeros(3*n,nk+4);
 t0 = ones(3*n,1);
 t1 = trm;
 t2 = trm.^2;
 spreg_m(:,1) = t0;
 spreg_m(:,2) = t1;
 spreg_m(:,3) = t2;
 i = 1; 
 while i <= nk+1;
	tmp = t1-nks(i);
	spreg_m(:,3+i) = (0.5*(tmp.^3)).*(t1 > nks(i));
    i=i+1; 
 end
 smt = spreg_m*b;
 %@ b';;sumc(llf); @
 %retp(sqt,smt);
 end            