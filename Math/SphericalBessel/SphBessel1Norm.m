% dimension order: n -> x
% =========================================================================
% 计算归一化第一类球Bessel函数
% -------------------------------------------------------------------------
% 输入
%	n		--- 阶数；必须在第1维
%	z		--- 宗量；必须在第2维及以上；绝对值必须小于1400
% =========================================================================
function jbar = SphericalBesselJNorm(n, z, varargin)
	p = inputParser;
	addParameter(p, 'isPrintInfo', 1);
	parse(p, varargin{:});
	ip = p.Results;

    % N_MAX = max(n(:)) + 2;
    z_row = z(:).';
    % 精度
    % prec = 1e-4;

    % N_MAX = max(ceil(1/2 * sqrt(max(abs(z_row))^2 * prec^(-1/3) + 9)), max(n(:))+3);
    N_MAX = max(n(:)) + round(2*abs(max(z(:))))+9;
    % N_MAX = max(n(:)) + 1e2;
    n_buf = (0:N_MAX).';
    % if n == 0
        % jbar = sinc(z/pi);
        % return
    % end

    % [~, idx_z_max] = max(abs(z_row));
    % z_max = z_row(idx_z_max);
    jbar = 0 * n_buf .* z_row + 1;         
    % T = 1e-30*z_max; % T=jbar(N+1)/jbar(N)
    % C = 1e-30*z_max; % used to compute T
    % D = 0*z_max; % used to compute T
    % k = N_MAX;
    % isExit = 0;
    % jbar(end,:) = 1e-300;
    % jbar(end-1,:) = 1e-300;

    % while isExit==0
        % D=(-4*(k+1)^2+1)./z.^2.* (1+D);    
        % for i=1:length(z)
            % if D(i)==0 || isnan(D(i)) || abs(D(i))==inf; D(i)=1e-30; end
        % end
        % D=1./D;
        % C=(-4*(k+1)^2+1)./z.^2.*(1+1./C);
        % for i=1:length(z)
            % if C(i)==0 || isnan(C(i)) || abs(C(i))==inf; C(i)=1e-30; end
        % end
        % Delta=C.*D;
        % T=T.*Delta;
        % if max(abs(1-Delta(:)))<1e-14 
            % isExit=1;
        % end
        % k=k+1;
    % end
    % arbitrary starting value to start downward recurrence:    
    % jbar(N_MAX+1, :) = ones(size(z_buf));
    % jbar(N_MAX, :) = ones(size(z_buf));
    % jbar(N_MAX,:) = jbar(N_MAX+1,:)./T;
        for nn = fliplr(1:N_MAX-1)
            jbar(nn, :) = jbar(nn+1, :) - z_row.^2.*jbar(nn+2,:)./(4*(nn+1)^2-1);
        end
    % jbar(1:N_MAX-1, :) = jbar(2:N_MAX, :) ...
    % 	- z_row.^2 ./ (4*((2:N_MAX).').^2 -1) .* jbar(3:N_MAX+1, :);
    %now all jbar satisfy the recurrence, just need to scale them to have the
    %correct values:
    % if abs(z)<=1400 % for Re{x}<=1400 (and Im{x}<700), can compute jbar for all n
    % if isnan(jbar(1))==false && jbar(1)~=0 && abs(jbar(1))~=inf 
        j0 = sinc(z_row/pi) ./ jbar(1,:);
        % for k=0:N_MAX
            % jbar(k+1,:)=jbar(k+1,:).*j0;
        % end
        jbar = jbar .* j0;
    % elseif N_MAX>abs(4/3*z+1) %for x>1400 have to discard jbar for n<~x because some values overflow 
        % - get the scale factor from the unnormalised Bessel function for n=4/3 x
        % n=floor(abs(4/3*z));
        % jn0=sqrt(pi/2/z)*besselj(n+1/2,z)*prod([3:2:2*n+1]/z) ./jbar(n+1,:);
        % for k=N_MAX:-1:n
            % jbar(k+1,:)=jbar(k+1,:).*jn0;
        % end    
    % end

    jbar = jbar(n+1,:);
    % jbar(:, abs(z_row)<1e-300) = 1;
    jbar = reshape(jbar, size(0*n.*z));

end
