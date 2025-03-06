% Intuition: illustration of the identification coming from non-gaussianity
% in a theoretical Demand-Supply example.
% This script produces Figure 2 and figure A.1 in the Online Appendix.
clear variables, close all

rotmat = @(a) [cos(a) sin(a); -sin(a) cos(a)];

% Theoretical example
rho = 0.2;
V = [1 rho; rho 1];
Ch = chol(V);
if rho>0
    alimits = [atan(Ch(1,2)/Ch(2,2)); pi/2];
    w1 = 0.9; w2 = 0.1;
    a1 = [w1 1-w1]*alimits;
    a2 = [w2 1-w2]*alimits;
end    
C1 = rotmat(a1)'*Ch;
C2 = rotmat(a2)'*Ch;

slopes1 = C1(:,2)./C1(:,1)
slopes2 = C2(:,2)./C2(:,1)


% plotting
maxabs = 6;
xmin = -maxabs; xmax = maxabs; ymin = -maxabs; ymax = maxabs;
pos = [5,3,20,8];

b_slopes = 1;
fh = figure('Units','centimeters','Position',pos);
subplot(1,2,1)
xline(0)
yline(0)
% add the slopes of the shocks
if b_slopes
    C = C1*100;
    line([-C(1,1) C(1,1)],[-C(1,2) C(1,2)],'LineStyle','-','Color','k','LineWidth',1.5)
    line([-C(2,1) C(2,1)],[-C(2,2) C(2,2)],'LineStyle',':','Color','k','LineWidth',1.5)
    text(0.75*xmin,C(1,2)/C(1,1)*0.75*xmin,'Demand','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','Rotation',-9)
    text(0.25*xmax,C(2,2)/C(2,1)*0.25*xmax,'Supply','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','Rotation',73)
end
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('\Delta Quantity')
ylabel('\Delta Price')
subplot(1,2,2)
xline(0)
yline(0)
% add the slopes of the shocks
if b_slopes
    C = C2*100;
    line([-C(1,1) C(1,1)],[-C(1,2) C(1,2)],'LineStyle','-','Color','k','LineWidth',1.5)
    line([-C(2,1) C(2,1)],[-C(2,2) C(2,2)],'LineStyle',':','Color','k','LineWidth',1.5)
    text(0.1*xmin,C(1,2)/C(1,1)*0.1*xmin,'Demand','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','Rotation',-82)
    text(0.75*xmax,C(2,2)/C(2,1)*0.75*xmax,'Supply','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','Rotation',20)
end
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('\Delta Quantity')
ylabel('\Delta Price')
exportgraphics(fh, 'models.pdf')

%% simulate data from Models 1 and 2
ndraws = 1000;
spec = 'pdmt';
specs = {'g','t','pdmt','mvt'};
rng(0)
Z = randn(ndraws,2);
for i = 1:length(specs)
spec = specs{i};
switch spec
    case 'g'
        U = Z;
        maxabs = 6;
    case 't'
        v = 1.5;
        q = chi2rnd(v,ndraws,2);
        U = Z.*(q.^-0.5).*(v^.5);
%        U = trnd(v,ndraws,2);
        U = U./std(U);
        maxabs = 6;
    case 'mvt'
        rng(12345)
        v = 1.5;
        q = chi2rnd(v, ndraws, 1);
        U = Z.*repmat((q.^-0.5).*(v^.5), 1, 2);
        U = U./std(U);
        maxabs = 6;
    case 'pdmt'
        rng(123468);
        v0 = .85;
        vn = 1.5-v0;
        q0 = chi2rnd(v0, ndraws, 1);
        qn = chi2rnd(vn, ndraws, 2);
        U = Z.*((q0 + qn).^-0.5).*sqrt(v0 + vn);
        U = U./std(U);
    case 'u'
        U = rand(ndraws, 2) - 0.5;
        U = U./std(U);
    case 'tr'
        pd = makedist('Triangular','a',-1,'b',0,'c',1);
        U = random(pd, ndraws, 2);
        U = U./std(U);
end

M1 = U*C1;
M2 = U*C2;

xmin = -maxabs; xmax = maxabs; ymin = -maxabs; ymax = maxabs;
for b_slopes = [0, 1]
fh = figure('Units','centimeters','Position',pos);

x = M1(:,1); y = M1(:,2);
subplot(1,2,1)
scatter(x,y)
xline(0)
yline(0)
% add the slopes of the shocks
if b_slopes
    C = C1*100;
    line([-C(1,1) C(1,1)],[-C(1,2) C(1,2)],'LineStyle','-','Color','k','LineWidth',1.5)
    line([-C(2,1) C(2,1)],[-C(2,2) C(2,2)],'LineStyle',':','Color','k','LineWidth',1.5)
end
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('\Delta Quantity')
ylabel('\Delta Price')

x = M2(:,1); y = M2(:,2);
subplot(1,2,2)
scatter(x,y)
xline(0)
yline(0)
% add the slopes of the shocks
if b_slopes
    C = C2*100;
    line([-C(1,1) C(1,1)],[-C(1,2) C(1,2)],'LineStyle','-','Color','k','LineWidth',1.5)
    line([-C(2,1) C(2,1)],[-C(2,2) C(2,2)],'LineStyle',':','Color','k','LineWidth',1.5)
end
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('\Delta Quantity')
ylabel('\Delta Price')

if strcmp(spec, 't')
    fname = [spec replace(num2str(v),'.','_')];
else
    fname = spec;
end
if b_slopes
    fname = [fname '-slopes'];
end
exportgraphics(fh, [fname '.pdf'], 'ContentType', 'vector')
end

% Report the likelihood across rotations
if strcmp(spec, 't')
    npoints = 50;
    aa = 0:pi/2/(npoints-1):pi/2;
    ll_t = nan(1,npoints);
    ll_g = nan(1,npoints);
    for i = 1:npoints
        Q = rotmat(aa(i));
        Wi = inv(Ch)*Q;
        ll_t(i) = loglik_MW_iidstudent(M1, Wi, 1.5);
        ll_g(i) = loglik_MW_normal(M1, Wi);
    end
    fh = figure('Units','centimeters','Position',[5,3,18,6]);
    subplot(1,2,1)
    hold on
    plot(aa,ll_t)
    xline(a1)
    xline(a2)
    text(a1, quantile(ylim,0.55), 'Model 1','HorizontalAlignment','center')
    text(a2, quantile(ylim,0.55), 'Model 2','HorizontalAlignment','center')
    ylabel('Log-Likelihood')
    xlabel('\alpha')
    xticks([0 pi/4 pi/2])
    xticklabels({'0','\pi/4','\pi/2'})
    title('Independent Student-t likelihood')
    subplot(1,2,2)
    hold on
    plot(aa,ll_g)
    ylim([-1000, -960])
    xline(a1)
    xline(a2)
    text(a1, quantile(ylim,0.55), 'Model 1','HorizontalAlignment','center')
    text(a2, quantile(ylim,0.55), 'Model 2','HorizontalAlignment','center')
    ylabel('Log-Likelihood')
    xlabel('\alpha')
    xticks([0 pi/4 pi/2])
    xticklabels({'0','\pi/4','\pi/2'})
    title('Gaussian likelihood')
    exportgraphics(fh, 'likelihoods.pdf')

    lr = -2*(loglik_MW_iidstudent(M1, inv(Ch)*rotmat(a2), 1.5) - loglik_MW_iidstudent(M1, inv(Ch)*rotmat(a1), 1.5))
    cutoff1pct = chi2inv(0.99,1)
end


end