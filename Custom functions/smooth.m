function y_smthd = smoothing(x, y, sp)

f = fit(x, y, 'smoothingspline', 'smoothingparam', sp);

y_smthd = f(x);

end