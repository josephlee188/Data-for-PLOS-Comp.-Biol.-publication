function value = dotp(vec1, vec2, dim)

% vec1 , vec2 must have same dimensions

sz = size(vec1);
% inds is a cell whose each element contains indecies (1~length) for the
% dimension.

inds = repmat({1},1,ndims(vec1)); % Initialize
for i=1:ndims(vec1)
    inds{i} = 1:sz(i);
end


sz0 = sz;
sz0(dim) = 1;
value = zeros(sz0);
for d = 1:sz(dim)
    inds{dim} = d;
    x1 = vec1(inds{:});
    x2 = vec2(inds{:});
    value = value + x1 .* x2;
end




% inds{dim} = 1;
% x1 = vec1(inds{:});
% x2 = vec2(inds{:});
% inds{dim} = 2;
% y1 = vec1(inds{:});
% y2 = vec2(inds{:});
% 
% value = x1.*x2 + y1.*y2;




end