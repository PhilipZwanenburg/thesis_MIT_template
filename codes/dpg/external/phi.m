function [Chi_v] = phi(P,XYZr)

% Returns the matrix 'Chi_v' representing the orthonormal basis functions 
% evaluated at the provided quadrature nodes for polynomial order P. Works
% for dimensions d = 1-3.

[Nn d] = size(XYZr);
N = P+1;
Nb = N^d;

Chi_v = zeros(Nn,Nb);

IndChi = 0;
if (d == 1)
    for i = 1:N
        IndChi = IndChi+1;
        Chi_v(:,IndChi) = JacobiP(XYZr(:,1),0,0,i-1);
    end
elseif (d == 2)
    for j = 1:N
    for i = 1:N
        IndChi = IndChi+1;
        Chi_v(:,IndChi) = JacobiP(XYZr(:,1),0,0,i-1).*...
                          JacobiP(XYZr(:,2),0,0,j-1);
    end
    end
elseif (d == 3)
    for k = 1:N
    for j = 1:N
    for i = 1:N
        IndChi = IndChi+1;
        Chi_v(:,IndChi) = JacobiP(XYZr(:,1),0,0,i-1).*...
                          JacobiP(XYZr(:,2),0,0,j-1).*...
                          JacobiP(XYZr(:,3),0,0,k-1);
    end
    end
    end
end

return;
