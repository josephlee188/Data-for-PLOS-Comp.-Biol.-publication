function [result, ctime] = cme2(position)

T = size(position, 1);
v = diff(position, 1);
V = sqrt(v(:,1).^2 + v(:,2).^2);
V_temp = V(~isnan(V));
el_temp = cumsum(V_temp);
el_temp = cat(1, 0, el_temp);
f = find(~isnan(V));
f = [f(1) ; (f + 1)];
el = nan(T, 1);
el(f) = el_temp;

result = nan(T, 1);
x = [];
for dt=0:T-1
    t = 1:T-dt;
    D2 = position(t+dt,:)-position(t,:);
    D2 = sqrt(squeeze(D2(:,1).^2+D2(:,2).^2));
    dL = el(t+dt)-el(t);
    ratio = D2 ./ dL;
    result(dt+1) = nanmean(ratio, 1);
    x(dt+1) = nanmean(dL);
end

% x = 0:T-1;
y = result(~isnan(result));
x = x(~isnan(result))';
% f = fittype('A*exp(-x*a) + B*exp(-x*b) + (1-A-B)', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'A', 'a','B', 'b'});
f = fittype('A*exp(-x*a) + B*exp(-x*b) + (1-A-B)', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'A', 'a', 'B', 'b'});
% myfit = fit(x,y, f , 'startpoint', [1, 1, 1, 1], 'lower', [0, 0.001, 0, 0.001], 'upper', [1, 20, 1, 20]);
myfit = fit(x,y, f , 'startpoint', [1, 1, 1, 1], 'lower', [0, 0.001, 0, 0.01], 'upper', [1, 100, 1, 100]);
plot(myfit, x, y); legend('off')
ctime = coeffvalues(myfit);



end