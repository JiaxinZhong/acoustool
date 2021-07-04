% 计算归一化球Bessel函数，见Zhong2020JASA式(26)
% dimorder: n -> x
function hbar = SphericalHankelNorm(n, z, varargin)
	p = inputParser;
	addParameter(p, 'isPrintInfo', 1);
	parse(p, varargin{:});
	ip = p.Results;
	z_row = z(:).';
	N_MAX = max(n)+1e2;

	% if abs(imag(z_row) < 1e-200)
	%     jbar = cal_sphBesselNorm(n,z_row);
	%     ybar = cal_sphNeumannNorm(n,z_row);
	%     hbar = ybar + 1i*pi*exp((2*n+1).*log(z_row/2) ...
	%         -gammaln(n+3/2) - gammaln(n+1/2)) .*jbar;
	%     return
	% end

	hbar = zeros(N_MAX+1,length(z_row));
	% if abs(x) <= 1400 % for x<1400 (and Im{x}<700), can compute hb via recurrence for all n 
		hbar(1,:) = exp(1i*z_row);    
		hbar(2,:) = 1i*(sin(z_row)-z_row.*cos(z_row)) + z_row.*sin(z_row)+cos(z_row);  
		for nn = 1:N_MAX-1
			hbar(nn+2,:) = hbar(nn+1,:) - z_row.^2./(4*nn^2-1).*hbar(nn,:);
		end
	%     T(2,:)=hb(2,:)./hb(1,:); % alternate compoutation using a continued fraction - gives similar results 
	%     for k=1:N-1
	%         T(k+2,:)=1-x.^2/(4*k^2-1)./T(k+1,:);
	%     end
	%     for k=2:N
	%         hb(k+1,:)=T(k+1,:).*hb(k,:);
	%     end
	% else % for x>1400 have to discard hb for n<~x because some values overflow - instead start the recurrence from n=4/3 x
	%     n0=floor(abs(4/3*x));
	%     hbar(n0+1,:)=sqrt(pi/2/x)*(besselj(n0+1/2,x) + 1i*bessely(n0+1/2,x)) *x/prod([1:2:2*n0-1]/x);
	%     hbar(n0+2,:)=sqrt(pi/2/x)*(besselj(n0+3/2,x) + 1i*bessely(n0+3/2,x)) *x/prod([1:2:2*n0+1]/x);
	%     for nn=n0+1:N-1
	%         hbar(nn+2,:)=hbar(nn+1,:) - x.^2./(4*nn^2-1).*hbar(nn+1,:);
	%     end
	% end
	% hbar(:, abs(z_row)<1e-200) = 0;
	% hbar(1, abs(z_row)<1e-200) =1;
	hbar = hbar(n+1,:);
	hbar = reshape(hbar, size(0*n.*z));

% 	if ip.isPrintInfo
% 		if check_infnan(hbar, 'mode', 'mute')
% 			warning('There exists inf or nan!\n');
% 		end
% 	end
end
