function ll = loglpdmt(Y, W, q0, qn, v0, vn)
% Evaluate the log likelihood of the PDMT model
% conditional on the q's

term1 = -0.5*numel(Y)*log(2*pi) + size(Y,1)*log(abs(det(W)));

Q = (qn + q0)./(vn + v0);

Z = Y*W.*sqrt(Q);

ll = term1 + 0.5*sum(log(Q),'all') - 0.5*sum(Z.^2,'all');