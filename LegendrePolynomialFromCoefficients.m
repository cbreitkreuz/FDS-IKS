function f = LegendrePolynomialFromCoefficients(coefficients,nlat)
% Written by C.Breitkreuz (last modified 31.01.2019)
% Computes LegendrePolynomial of degree N, from coefficients

% Input:
% coefficients - coefficients of Legendre polynomial
% nlat - number of latitude coordinates

% Output:
% f - Sum of Legendre polynomials

N = length(coefficients);

x = linspace(-1, 1, nlat)' ;


% Legendre Polynomial
% P(n) = P_n-1, e.g. P(1) = P0

P = zeros(nlat,N);

P(:,1) = ones(1,nlat);

if N > 1
    P(:,2) = x;
    
    if N > 2
        % Compute Legendre plynomial with recursion formel (see Wikipedia)
        % Careful, it's extra confusing because P(:,n) = P_n-1
        for n = 2:N-1
            P(:,n+1) = ((2*(n-1) + 1) .* x .* P(:,n) - (n-1) .* P(:,n-1))./(n-1+1);
        end
    end
end

% If N = 3 this is equivalent to:
% P0 = ones(1,length(x))';
% P1 = x;
% P2 = 0.5.*(3.*x.^2 -1);
% P3 = 0.5.*(5.*x.^3 -3.*x);

f = zeros(nlat,1);

for lat = 1:nlat
    f(lat) = coefficients' * P(lat,:)';
end

end

