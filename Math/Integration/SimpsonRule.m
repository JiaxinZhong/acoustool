function res = SimpsonRule(integrand,  stepsize, dim)
sz = size(integrand);

inds = repmat({1}, 1, ndims(integrand));

for i = 1:length(sz)
    inds{i} = 1:sz(i);
end


inds1 = inds;
inds1{dim} = 1;

% inds0 = inds;
% inds0{dim} = 2;
% h = integrand(inds0{:}) - integrand(inds1{:});

inds2 = inds;
inds2{dim} = 2:2:sz(dim)-1;
inds3 = inds;
inds3{dim} = 3:2:sz(dim)-2;
inds4 = inds;
inds4{dim} = sz(dim);

res = stepsize/3.*(integrand(inds1{:}) + 4 * sum(integrand(inds2{:}), dim) ...
    + 2 * sum(integrand(inds3{:}), dim) ...
    + integrand(inds4{:}));