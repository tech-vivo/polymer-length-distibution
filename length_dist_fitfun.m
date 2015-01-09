% Will McFadden (wmcfadden)
% fitting function referenced in pl_fitting

function ret = length_dist_fitfun(q, dat)
    cc = jet(101);
    nmax = 100;
    t = dat{1};
    
    
    nuc = q(1);
    chnk = 0;
    for i=1:nmax-1
        chnk = [chnk; chnk(end)+50*i];
    end
    chnk = chnk + round(nuc);
    kon = q(2);
    koff = q(3);
    knuc = q(4)*10^-9;
    ksev = q(5)*10^-8;
    c = q(6)*100;
    kph = q(7)*10^-5;
    
    
    fprintf('%d %d %d %d %d %d %d \n', q(1), q(2), q(3), q(4), q(5), q(6), q(7));

    dT = 10;

    ret = [];
    for act_0 = [1 2 3 4 5]
        tic %odeset('NonNegative',1:nmax)
        [t, a] = ode15s(@length_dist_ode,t,zeros(nmax,1),odeset('NonNegative',1:nmax,'AbsTol',1e-16),act_0,chnk,kon,koff,knuc,nuc,ksev);
        toc
%         figure;
%         hold off
%         plot(chnk/370, a(end,:));
%         hold on
        asum = zeros(size(a,1),1);
        bsum = zeros(size(a,1),1);
        cc = jet(size(a,1));
        for ind = 2:size(a,1)
%             if(mod(ind-1,dT/10)==0)
                a2 = a(ind,:)';
                aa = nuc*a2(1)+sum(diff(chnk).*chnk(1:end-1).*a2(1:end-1) + chnk(1:end-1).*diff(a2).*(diff(chnk)+1)/2 + a2(1:end-1).*diff(chnk).*(diff(chnk)+1)/2 + diff(a2).*(diff(chnk)+1).*(2*diff(chnk)+1)/6);
                bsum(ind) = asum(ind-1)+(bsum(ind-1)-asum(ind-1))*exp(-kph*10);
                asum(ind) = aa;
%                 plot(chnk/370, a2,'.','color',cc(ind,:));
%             end
        end
        ret = [ret c*asum.*exp(-kph*t)];
    end
end