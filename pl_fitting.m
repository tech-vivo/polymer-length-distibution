% Will McFadden (wmcfadden)
% fits simulation results to actin polymerization data
% given input of time t and fluorescence R;
% assumes it is given actin at concentrations 1, 2, 3, 4, and 5 uM

fitval = R;
fitdat = {t};
param =[3.0590    5.1676    0.1193    6.5802    0.1181    4.9451    4.7469];
param = [3.006247e+00 5*2.120564e+00 5*7.507815e-02 (1/5)*1.534126e+01 2.762542e-01 5.020619e+00 4.957969e+00 ];
% c0 = 5;
% kph0 = 5;
% 
% nuc0 = 3;
% kon0 = 10;
% koff0 = 1;
% knuc0 = 4;
% ksev = 1;

rjump = 0.1;
sto = [];
for i = 1:1
    
%     param = [nuc0,kon0,koff0,knuc0,ksev,c0,kph0];
    lb = [2.5, 0, 0, 0, 0, 2, 0];
    up = [4.1, 16, 4, 1000, 10, 7, 100];
    [sol,MSE,residual,exitflag,output,lambda,J] = lsqcurvefit(@length_dist_fitfun,param,fitdat,fitval,lb,up,optimset('FinDiffRelStep',0.001));

    
    nuc = sol(1);
    kon = sol(2);
    koff = sol(3);
    knuc = sol(4);
    ksev = sol(5);
    c = sol(6);
    kph = sol(7);
    
    
    sto = [sto; nuc kon koff knuc c kph];
    c0 = c*(1 + rjump*(rand-0.5));
    kph0 = kph*(1 + rjump*(rand-0.5));
    nuc0 = nuc*(1 + rjump*(rand-0.5));
    kon0 = kon*(1 + rjump*(rand-0.5));
    koff0 = koff*(1 + rjump*(rand-0.5));
    knuc0 = knuc*(1 + rjump*(rand-0.5));
    rjump = rjump*0.5;
end
% subplot(2,2,1)
% plot(sto(:,1))
% ylabel('k_{on}')
% subplot(2,2,2)
% plot(sto(:,2))
% ylabel('k_{off}')
% subplot(2,2,3)
% plot(sto(:,3))
% ylabel('k_{ph}')
% subplot(2,2,4)
figure;
plot(t,R,'.')
hold on
plot(t,length_dist_fitfun(sol,fitdat),'r')
ylabel('final fit')
finalanswer = full([c; kph; nuc; kon; koff; knuc])