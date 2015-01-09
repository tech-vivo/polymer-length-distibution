for i=1:size(R,2)
    R(:,i)=R(:,i)-R(1,i);
end
fitval = R;
fitdat = {t};
% param = [2.9461    4.4239    0.1688    7.3842    0.1969    5.1176    4.9803];
param = [3.006247e+00 5*2.120564e+00 5*7.507815e-02 (1/5)*1.534126e+01 2.762542e-01 5.020619e+00 4.957969e+00 ];

y = length_dist_fitfun(param,fitdat);

figure;
plot(t,R,'.')
hold on
plot(t,y,'r')
ylabel('final fit')
