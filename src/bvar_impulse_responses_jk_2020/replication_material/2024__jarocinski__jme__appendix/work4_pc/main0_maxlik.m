% Estimate Y = UC, Student-t, unknown v
% where
% Y - observable data, i.i.d.
% U - i.i.d. Student-t WITH UNKNOWN v
% C - matrix with the impacts
% We will parameterize in terms of W and z, such that
% YW = U, then recover C = inv(W)
% z = log(v), then recover v = exp(z)
%
% Maximum Likelihood
% Depends: nlogl_iidstudent_Wzn.m (->kr.m), normalize1.m
clear all, close all

% read data
%data_pcx3
%data_pcx3_SP500
%data_pcxprime4
%data_pcxprime3
data_pcxprime5
%data_pcxbis5
%data_pcxbis6

[T, N] = size(Y);
nv = N; % 1 or N

out_path = 'out/';
mkdir(out_path)

figure()
for n = 1:N
subplot(ceil(sqrt(N)),floor(sqrt(N)),n)
plot(tab.date, Y(:,n))
title("Y" + n)
end

% Maximum likelihood / posterior mode estimation of W(->C), z(->v)
options = optimoptions('fminunc');
options.Display = 'iter';
options.MaxFunctionEvaluations = 1e4*N^2;
options.OptimalityTolerance = 1e-8;
options.Algorithm = 'trust-region';
options.SpecifyObjectiveGradient = true;
% starting point
W = inv(chol(cov(Y)) + randn(N,N));
par0 = [W(:); repmat(log(3),N,1)];

% maximum likelihood
fun = @(par) nlogl_iidstudent_Wzn(Y, reshape(par(1:N^2),N,N), par(N^2+1:end));
[parmaxlik,fval,exitflag,output,grad,hessian] = fminunc(fun, par0, options);
W = reshape(parmaxlik(1:N^2), N, N);
C = inv(W);
v = exp(parmaxlik(N^2+1:end));

% Normalize
P = eye(N); P = diag(sign(sorder))*P(abs(sorder),:);

%P*C
Wnormalized = W*P';
vnormalized = abs(P)*v;


% Minimize again to get the reordered hessian
par0 = [Wnormalized(:); log(vnormalized)];
[parmaxlik,fval,exitflag,output,grad,hessian] = fminunc(fun, par0, options);

maxl.Y = Y;
maxl.ll = -fval;
maxl.parmaxlik = parmaxlik;
maxl.hessian = hessian;
maxl.v = exp(parmaxlik(N^2+1:end));
maxl.W = reshape(parmaxlik(1:N^2), N, N);
maxl.C = inv(maxl.W);
maxl.U = Y*maxl.W;
maxl.U1s = maxl.U./std(maxl.U);
maxl.C1s = diag(std(maxl.U))*maxl.C;
maxl.vdec = (maxl.C1s.^2)./sum(maxl.C1s.^2);
if N==4, temp = [1 2 3 2]; else, temp = 1:N; end
temp = diag(diag(maxl.C(:,temp)));
maxl.C1bp = temp\maxl.C;
maxl.U1bp = maxl.U*temp;

disp('C'), disp(maxl.C)
disp('C1s'), disp(maxl.C1s)
disp('C1bp'), disp(maxl.C1bp)
disp('vdec'), disp(maxl.vdec)
disp('v'), disp(maxl.v)

% save data, maxl, shocks
save(out_path + "data", "tab");
save(out_path + "maxl", "maxl");
u = round(maxl.U, 5);
writetable([tab(:,'date'), array2table(u)], out_path + "U.csv");
u = round(maxl.U1s, 5);
Utable = [tab(:,'date'), array2table(u)];
writetable(Utable, out_path + "U1s.csv");
u = round(maxl.U1bp, 5);
Utable = [tab(:,'date'), array2table(u)];
writetable(Utable, out_path + "U1bp.csv");

% Variance of C, v
JacobianWC = -kron(maxl.W',maxl.W);
Jacobianzv = diag(maxl.v.^-1);
JacobianWzCv = blkdiag(JacobianWC, Jacobianzv);
asyvarCv = inv(JacobianWzCv'*hessian*JacobianWzCv);
asystdCv = sqrt(diag(asyvarCv));
Cstd = reshape(asystdCv(1:N^2),N,N);
vstd = asystdCv(N^2+1:end);
disp(array2table([maxl.v vstd], ...
    "VariableNames", ["MLE","asy-std"], "RowNames", "v"+(1:N)))

% Plot the response of the Y variables (not very interpretable...)
Cm = maxl.C;
Cl = Cm - 2*Cstd;
Cu = Cm + 2*Cstd;
fh = plot_resp(diag(std(maxl.U))*Cm, diag(std(maxl.U))*Cl, diag(std(maxl.U))*Cu, ymaturities, ynames);
exportgraphics(fh, out_path + "C1smaxlik_band.pdf")

% Plot the response of the XZ variables
CCm = Cm*Y2XZ;
CCl = Cl*Y2XZ;
CCu = Cu*Y2XZ;

fh = plot_resp(diag(std(maxl.U))*CCm, diag(std(maxl.U))*CCl, diag(std(maxl.U))*CCu, xzmaturities, xznames);
exportgraphics(fh, sprintf('%sCC1smaxlik_band.pdf', out_path))
