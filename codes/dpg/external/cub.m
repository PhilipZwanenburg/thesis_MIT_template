function [xir,W,Nn] = cub(P,d,type)
% P = 2; d = 1; type = 'ES';

% Return nodes/weights for tensor product quadrature depending on 
% quadrature type.

N = P+1;

if (strcmp(type,'GL'))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function determines the abscisas (x) and weights (w)  for the      %
% Gauss-Legendre quadrature, of order n>1, on the interval [-1, +1].      %
%   Unlike many publicly available functions, 'GaussLegendre_2' is valid  %
%   for n>=46. This is due to the fact that 'GaussLegendre_2' does not    %
%   rely on the build-in Matlab routine 'roots' to determine the roots of %
%   the Legendre polynomial, but finds the roots by looking for the       %
%   eigenvalues of an alternative version of the companion matrix of the  %
%   n'th degree Legendre polynomial. The companion matrix is constructed  %
%   as a symmetrical matrix, guaranteeing that all the eigenvalues        %
%   (roots) will be real. On the contrary, the 'roots' function uses a    %
%   general form for the companion matrix, which becomes unstable at      %
%   higher orders n, leading to complex roots.                            %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Geert Van Damme
% geert@vandamme-iliano.be
% February 21, 2010    
% http://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-
% laguerre-and-hermite-gauss-quadrature/content/GaussLegendre.m

% Building the companion matrix CM
    % CM is such that det(xI-CM)=P_n(x), with P_n the Legendre polynomial
    % under consideration. Moreover, CM will be constructed in such a way
    % that it is symmetrical.
    i   = 1:P;
    a   = i./sqrt(4*i.^2-1);
    CM  = diag(a,1) + diag(a,-1);

% Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=P_n(x), the abscissas are the roots of the
    %   characteristic polynomial, i.d. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
    [V L]   = eig(CM);
    [x ind] = sort(diag(L));
    V       = V(:,ind)';
    w       = 2 * V(:,1).^2;
elseif (strcmp(type,'GLL'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% lglnodes.m                                                              %
%                                                                         %
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL          %
% Vandermonde matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x).     %
% Useful for numerical integration and spectral methods.                  %
%                                                                         %
% Reference on LGL nodes and weights:                                     %
% C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods %
% in Fluid Dynamics," Section 2.3. Springer-Verlag 1987                   %
%                                                                         %
% Written by Greg von Winckel - 04/17/2004                                %
% Contact: gregvw@chtm.unm.edu                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Truncation + 1
    P1=P+1;

    % Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    if (P == 0); 
        disp('Cannot use GLL nodes of order P0');
        pause
    else
        x=cos(pi*(0:P)/P)';
    end

    % The Legendre Vandermonde Matrix
    V=zeros(P1,P1);

    % Compute P_(N) using the recursion relation
    % Compute its first and second derivatives and 
    % update x using the Newton-Raphson method.

    xold=2;

    while max(abs(x-xold))>eps
        xold = x;
        V(:,1)=1;    V(:,2)=x;

        for k=2:P
            V(:,k+1)=( (2*k-1)*x.*V(:,k)-(k-1)*V(:,k-1) )/k;
        end
        x = xold-( x.*V(:,P1)-V(:,P) )./( P1*V(:,P1) );
    end

    w=2./(P*P1*V(:,P1).^2);
    
    x = (fliplr(x'))';
elseif (strcmp(type,'ES'))
    if (P == 0); x = 0;
    else         x = (-1:2/P:1)';
    end
    w = zeros(P+1,1);
end

xir = zeros(N^d,d);
W = zeros(N^d,1);

if (d == 1)
    xir = x;
    W = w;
elseif (d == 2)
    for i = 1:N
        xir((i-1)*N+1:i*N,:) = [x(:,1) x(i,1)*ones(N,1)];
        W((i-1)*N+1:i*N,:) = w(:,1)*w(i,1);
    end
elseif (d == 3)
    for j = 1:N
    for i = 1:N
        xir((j-1)*N^2+(i-1)*N+1:(j-1)*N^2+i*N,:) = ...
            [x(:,1) x(i,1)*ones(N,1) x(j,1)*ones(N,1)];
        W((j-1)*N^2+(i-1)*N+1:(j-1)*N^2+i*N,:) = ...
            w(:,1)*w(i,1)*w(j,1);
    end
    end
end

Nn = N^d;

if (strcmp(type,'ES'))
    % Find Connectivity for Tecplot
    if     (d == 1); NvnGs = 2;
    elseif (d == 2); NvnGs = 4;
    elseif (d == 3); NvnGs = 8;
    end

    Connect = zeros(P^d,NvnGs);
    
    IndC = 1;
    for i = 1:P
        nLINE = [i i+1];
        for j = 1:max(P*min(d-1,1),1)
            nQUAD = [nLINE+N*(j-1) fliplr(nLINE+N*j)];
            for k = 1:max(P*min((d-2),1),1)
                nHEX = [nQUAD+N^2*(k-1) nQUAD+N^2*k];
                
                Connect(IndC,:) = nHEX(1,1:NvnGs);
                IndC = IndC + 1;
            end
        end
    end
    
    W = Connect;
    Nn = [size(Connect,1) Nn];
end

return;
