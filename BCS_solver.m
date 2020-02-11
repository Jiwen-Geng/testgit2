function [weights,used,sigma2,errbars] = BCS_solver(PHI,t,sigma2,eta,TotalCnts)
% Name: BCS_solver. 
% Դ�� Shihao Ji, ECE, Duke University��BCS_fast_rvm�㷨���Ըó�����о���ȥ���������ҵ�һЩ���֡�
% Jiwen Geng;     18 July, 2019
% ������⣺  t = PHI*W + noise
%>>>>>>>����:
%   PHI: ͶӰ����
%   t: ѹ����֪��������
%  sigma2: ��ʼ����������
%       ����������ڻ���w������ϡ�裬���Ƽ� sigma2 = std(x)^2/1e2;
%       ���������w����ϡ�裬���Ƽ� sigma2 = std(x)^2/1e6
%  eta: �㷨ֹͣ���ޣ��Ƽ�ֵ��1e-8
%>>>>>>>�����
%  weights: ��Ȩϵ��
%  used: ��Ȩϵ����λ������
%  sigma2: ���¹��Ƶ���������
%  errbars: ��Ȩϵ���ı�׼ƫ��
% Updates: 
% 2019/8/9, �����������TotalCnts,���Կ���������ѭ��������Ĭ�� 1000��
if nargin <5
    TotalCnts = 1000;    % Ĭ��ѭ������
end

% find initial alpha
[N,M] = size(PHI);
PHIt = PHI'*t;
PHI2 = sum(PHI.^2)';
ratio = (PHIt.^2)./PHI2;
[maxr,index] = max(ratio);
alpha = PHI2(index)/(maxr-sigma2);
% compute initial mu, Sig, S, Q
phi = PHI(:,index);
Hessian = alpha + phi'*phi/sigma2;
Sig = 1/Hessian;
mu = Sig*PHIt(index)/sigma2;
left = PHI'*phi/sigma2;
S = PHI2/sigma2-Sig*left.^2;
Q = PHIt/sigma2-Sig*PHIt(index)/sigma2*left;
%
% ��������
for count = 1:TotalCnts
    s = S; q = Q;
    s(index) = alpha.*S(index)./(alpha-S(index));
    q(index) = alpha.*Q(index)./(alpha-S(index));
    theta = q.^2-s;
    
    % choice the next alpha that maximizes marginal likelihood
    ml = -inf*ones(1,M);
    ig0 = find(theta>0);
    % index for re-estimate
    [ire,foo,which] = intersect(ig0,index);
    if ~isempty(ire)
        Alpha = s(ire).^2./theta(ire);
        delta = (alpha(which)-Alpha)./(Alpha.*alpha(which));
        ml(ire) = Q(ire).^2.*delta./(S(ire).*delta+1)-log(1+S(ire).*delta);
    end
    % index for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
        ml(iad) = (Q(iad).^2-S(iad))./S(iad)+log(S(iad)./(Q(iad).^2));
    end
    is0 = setdiff([1:M],ig0);
    % index for deleting
    [ide,foo,which] = intersect(is0,index);
    if ~isempty(ide)
        ml(ide) = Q(ide).^2./(S(ide)-alpha(which))-log(1-S(ide)./alpha(which));
    end

    [ML(count),idx] = max(ml);
    % check if terminates?
    if count > 2 & abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta
        break;
    end

    % update alphas
    which = find(index==idx);
    if theta(idx) > 0
        if ~isempty(which) % re-estimate
            Alpha = s(idx)^2/theta(idx);
            Sigii = Sig(which,which); mui = mu(which); Sigi = Sig(:,which);
            delta = Alpha-alpha(which);
            ki = delta/(1+Sigii*delta);
            mu = mu-ki*mui*Sigi;
            Sig = Sig-ki*Sigi*Sigi';
            comm = PHI'*(phi*Sigi)/sigma2;
            S = S + ki*comm.^2;
            Q = Q + ki*mui*comm;
            %
            alpha(which) = Alpha;
        else % adding
            Alpha = s(idx)^2/theta(idx);
            phii = PHI(:,idx); Sigii = 1/(Alpha+S(idx)); mui = Sigii*Q(idx);
            comm1 = Sig*(phi'*phii)/sigma2;
            ei = phii-phi*comm1;
            off = -Sigii*comm1;
            Sig = [Sig+Sigii*comm1*comm1', off; off', Sigii];
            mu = [mu-mui*comm1; mui];
            comm2 = PHI'*ei/sigma2;
            S = S - Sigii*comm2.^2;
            Q = Q - mui*comm2;
            %
            index = [index;idx];
            alpha = [alpha;Alpha];
            phi = [phi,phii];
        end
    else
        if ~isempty(which) % deleting
            Sigii = Sig(which,which); mui = mu(which); Sigi = Sig(:,which);
            Sig = Sig-Sigi*Sigi'/Sigii; Sig(:,which) = []; Sig(which,:) = [];
            mu  = mu-mui/Sigii*Sigi; mu(which) = [];
            comm = PHI'*(phi*Sigi)/sigma2;
            S = S + comm.^2/Sigii;
            Q = Q + mui/Sigii*comm;
            %
            index(which) = [];
            alpha(which) = [];
            phi(:,which) = [];
        end
    end
end
weights	= mu;
used = index;
% re-estimated sigma2
sigma2 = sum((t-phi*mu).^2)/(N-length(index)+alpha'*diag(Sig)); 
errbars = sqrt(diag(Sig));
end