% Find posterior mode
% Run after main0_maxlik
% Depends: nlogpost_iidstudent_Wzn.m (->kr.m)

% prior for W
prior.uW = eye(N);
prior.Ginv = zeros(N^2);
% prior for v
prior.va = 2*ones(1,nv);
prior.vb = 5*ones(1,nv);
save(out_path + "prior", "prior")

% use the maximum likelihood estimate as the starting point
par = parmaxlik;

% find posterior mode
fun = @(par) nlogpost_iidstudent_Wzn(Y, reshape(par(1:N^2),N,N), par(N^2+1:end), prior);
[parmaxlik,fval,exitflag,output,grad,hessian] = fminunc(fun, par0, options);

W = reshape(parmaxlik(1:N^2), N, N);
C = inv(W);
v = exp(parmaxlik(N^2+1:end));

maxp.Y = Y;
maxp.ll = -fval;
maxp.parmaxlik = parmaxlik;
maxp.hessian = hessian;
maxp.v = exp(parmaxlik(N^2+1:end));
maxp.W = reshape(parmaxlik(1:N^2), N, N);
maxp.C = inv(maxp.W);
maxp.U = Y*maxp.W;
maxp.U1s = maxp.U./std(maxp.U);
maxp.C1s = diag(std(maxp.U))*maxp.C;
maxp.vdec = (maxp.C1s.^2)./sum(maxp.C1s.^2);
if N==4, temp = [1 2 3 2]; else, temp = 1:N; end
temp = diag(diag(maxp.C(:,temp)));
maxp.C1bp = temp\maxp.C;
maxp.U1bp = maxp.U*temp;

disp('C'), disp(maxp.C)
disp('C1s'), disp(maxp.C1s)
disp('C1bp'), disp(maxp.C1bp)
disp('vdec'), disp(maxp.vdec)
disp('v'), disp(maxp.v)

save([out_path 'maxp'],'maxp');



