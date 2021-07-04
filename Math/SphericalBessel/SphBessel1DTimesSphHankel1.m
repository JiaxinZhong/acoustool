%==========================================================================
% The derivative of the spherical Bessel function of the first kind times 
%   the spherical Hankel function of the first kind.
%--------------------------------------------------------------------------
% INPUT
%   N       --- the maximum degree
%   zj      --- the argument of the spherical Bessel function
%   zh      --- the argument of the spherical Hankel function
%   dimension: n => zj.*zh
%--------------------------------------------------------------------------
% OUTPUT
%   res     --- The result, where the first dimension is the results 
%               corresponding to the degree 0, 1, 2, ..., N
%==========================================================================

function res = SphBessel1DTimesSphHankel1(N, zj, zh, varargin)
    p = inputParser;
    addParameter(p, 'method', 'robust');
    parse(p, varargin{:});
    ip = p.Results;

    zj_row = zj + zj .* zh * 0;
    zj_row = zj_row(:).';
    zh_row = zh + zj .* zh * 0;
    zh_row = zh_row(:).';
    
    switch ip.method
        case 'direct'
            n = (0:N).';
            jD = SphericalBesselJD(n, zj_row);
            h = SphericalHankel(n, zh_row);
            res = jD .* h;
        case 'norm'
            n = (0:N+1).';
            h_bar = SphericalHankelNorm(n, zh_row);
            j_bar = SphericalBesselJNorm(n, zj_row);

            n = (0:N).';
            res = exp((n-1).*log(zj_row./zh_row) ...
                + log(h_bar(1:end-1) .* (n .* j_bar(1:end-1) ...
                - zj_row.^2 ./ (2*n+3) .* j_bar(2:end)))) ...
                ./ (1i * (2*n+1) .* zh_row.^2);
        case 'robust'
            n = (0:N+1).';
            h_bar = SphHankel1Norm(n, zh_row);
            j_bar = SphBessel1Norm(n, zj_row);
            h0 = SphHankel1(n, zh_row);
            j0 = SphBessel1(n, zj_row);

            n = (0:N).';
            res = exp((n-1).*log(zj_row./zh_row) ...
                + log(h_bar(1:end-1,:) .* (n .* j_bar(1:end-1,:) ...
                - zj_row.^2 ./ (2*n+3) .* j_bar(2:end,:)))) ...
                ./ (1i * (2*n+1) .* zh_row.^2);
            res0 = (n ./ zj_row .* j0(1:end-1,:) - j0(2:end,:)) .* h0(1:end-1,:);
            
            h_bar = h_bar(1:end-1,:);
            j_bar = j_bar(1:end-1,:);
            res(isnan(h_bar) | isinf(h_bar) | isnan(j_bar) | isinf(j_bar)) ...
                = res0(isnan(h_bar) | isinf(h_bar) | isnan(j_bar) | isinf(j_bar));
            res(isnan(res)) = 0;
            
            res(:, abs(zj_row)>1200 | abs(zh_row)>1200) ...
                = res0(:, abs(zj_row)>1200 | abs(zh_row)>1200);
    end
    
    res = reshape(res, size(n .* zj .* zh));
end
