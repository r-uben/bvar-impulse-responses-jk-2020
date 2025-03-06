function ll = loglik_MW_iidstudent(M, W, v)
% Evaluate log p(M|W,v)
% Model M = U*inv(W), U~i.i.d. student
U = M*W;
logpU = -(v+1)/2*log(U.^2/v + 1);
ll = size(M,1)*log(det(W)) + sum(sum(logpU));