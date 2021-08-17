function value = absvec(vector, dim)

sz = size(vector);
% inds is a cell whose each element contains indecies (1~length) for the
% dimension.
inds = repmat({1},1,ndims(vector)); % Initialize
for i=1:ndims(vector)
    inds{i} = 1:sz(i);
end

sz0 = sz;
sz0(dim) = 1;
value = zeros(sz0);
for d = 1:sz(dim)
    inds{dim} = d;
    x = vector(inds{:});
    value = value + x.^2;
end

value = sqrt(value);

% inds{dim} = 1;
% x = vector(inds{:});
% inds{dim} = 2;
% y = vector(inds{:});
% value = sqrt(x.^2 + y.^2);





% A = rand(100,100,10,10);
% dim = 4;
% sz = size(A);
% inds = repmat({1},1,ndims(A));
% inds{dim} = 1:sz(dim);
% A(inds{:})

end