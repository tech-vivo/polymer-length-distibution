% Will McFadden (wmcfadden)

function da = length_dist_ode(t,a,act_0,chnk,kon,koff,knuc,nuc,ksev)
    nmax = length(a);
    chall = [a(1); diff(chnk).*a(2:end)-diff(a).*(3*diff(chnk)+1)/2];
    asum = cumsum(chall);
    act = act_0-nuc*a(1)-sum(diff(chnk).*chnk(1:end-1).*a(1:end-1) + chnk(1:end-1).*diff(a).*(diff(chnk)+1)/2 + a(1:end-1).*diff(chnk).*(diff(chnk)+1)/2 + diff(a).*(diff(chnk)+1).*(2*diff(chnk)+1)/6);
    da = zeros(1,nmax)';
    
    da(1) = knuc*act^nuc - (koff+kon*act)*a(1) + koff*interp1q(chnk,a,chnk(1)+1)+2*ksev*(asum(end)-asum(1));
    for i = 2:nmax-1
        da(i) = kon*act*interp1q(chnk,a,chnk(i)-1) + koff*interp1q(chnk,a,chnk(i)+1)-(koff+kon*act+(chnk(i)-1)*ksev)*a(i) +2*ksev*(asum(end)-asum(i));
    end
    da(end) = kon*act*interp1q(chnk,a,chnk(i)-1) -(koff+kon*act+(chnk(i)-1)*ksev)*a(i); 
end