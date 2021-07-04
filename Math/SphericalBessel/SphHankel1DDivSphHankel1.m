%==========================================================================
% The derivative of the spherical Hankel function of the first kind divided 
%   by a spherical Hankel function of the first kind.
%--------------------------------------------------------------------------
% INPUT
%   N       --- the maximum degree
%   zd      --- the argument of the derivative of the spherical Hankel 
%               function
%   z0      --- the argument of the spherical Hankel function
%   dimension: n => zd .* z0
%--------------------------------------------------------------------------
% OUTPUT
%   res     --- The result, where the first dimension is the results 
%               corresponding to the degree 0, 1, 2, ..., N
%==========================================================================

function res = SphHankel1DDivSphHankel1(N, zd, z0, varargin)
    p = inputParser;
    addParameter(p, 'method', 'robust');
    parse(p, varargin{:});
    ip = p.Results;

    zd_row = zd + 0*zd.*z0;
    z0_row = z0 + 0*zd.*z0;
    zd_row = zd_row(:).';
    z0_row = z0_row(:).';
    
    switch ip.method
        case 'direct'
            n = (0:N).';
            hD = SphericalHankelD(n, zd_row);
            h = SphericalHankel(n, z0_row);
            res = hD ./ h;
        case 'norm'
            n = (0:N+1).';
            h_bar_zd = SphericalHankelNorm(n, zd_row);
            h_bar_z0 = SphericalHankelNorm(n, z0_row);

            n = (0:N).';
            res = exp((n+1).*log(z0_row./zd_row) ...
                + log(n .* h_bar_zd(1:end-1,:) ...
                - (2.*n+1) .* h_bar_zd(2:end,:)) ...
                - log( zd_row .* h_bar_z0(1:end-1,:)));

        case 'robust'
            n = (0:N+1).';
            h_bar_zd = SphHankel1Norm(n, zd_row);
            h_bar_z0 = SphHankel1Norm(n, z0_row);
            h_zd = SphHankel1(n, zd_row);
            h_z0 = SphHankel1(n, z0_row);

            n = (0:N).';
            res = exp((n+1).*log(z0_row./zd_row) ...
                + log(n .* h_bar_zd(1:end-1,:) ...
                - (2.*n+1) .* h_bar_zd(2:end,:)) ...
                - log( zd_row .* h_bar_z0(1:end-1,:)));
            res0 = (n ./ zd_row .* h_zd(1:end-1,:) - h_zd(2:end,:)) ...
                ./ h_z0(1:end-1,:);
            h_bar_z0 = h_bar_z0(1:end-1,:);
            h_bar_zd = h_bar_zd(1:end-1,:);
            res(isnan(h_bar_z0) | isinf(h_bar_z0) | isnan(h_bar_zd) | isinf(h_bar_z0)) ...
                = res0(isnan(h_bar_z0) | isinf(h_bar_z0) | isnan(h_bar_zd) | isinf(h_bar_z0));

            res(:, abs(zd_row)>1200 | abs(z0_row) > 1200) =...
                res0(:, abs(zd_row)>1200 | abs(z0_row) > 1200);
    end
    res = reshape(res, size(n .* zd .* z0));
end
