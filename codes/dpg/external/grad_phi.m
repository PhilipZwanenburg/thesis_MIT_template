function [GradChi_v] = grad_phi(P,XYZr)

% Returns the matrix 'GradChi_v' representing the gradient of the 
% orthonormal basis functions evaluated at the provided quadrature nodes 
% for polynomial order P. Works for dimensions d = 1-3.

[Nn d] = size(XYZr);
N = P+1;
Nb = N^d;

GradChi_v = zeros(Nn,Nb,d);

IndChi = 0;
if (d == 1)
    for i = 0:P
        IndChi = IndChi+1;
        GradChi_v(:,IndChi,1) = GradJacobiP(XYZr(:,1),0,0,i);
    end
elseif (d == 2)
    for j = 0:P
    for i = 0:P
        IndChi = IndChi+1;
        GradChi_v(:,IndChi,1) = GradJacobiP(XYZr(:,1),0,0,i).*...
                                    JacobiP(XYZr(:,2),0,0,j);
        GradChi_v(:,IndChi,2) =     JacobiP(XYZr(:,1),0,0,i).*...
                                GradJacobiP(XYZr(:,2),0,0,j);
    end
    end
elseif (d == 3)
    for k = 0:P
    for j = 0:P
    for i = 0:P
        IndChi = IndChi+1;
        GradChi_v(:,IndChi,1) = GradJacobiP(XYZr(:,1),0,0,i).*...
                                    JacobiP(XYZr(:,2),0,0,j).*...
                                    JacobiP(XYZr(:,3),0,0,k);
        GradChi_v(:,IndChi,2) =     JacobiP(XYZr(:,1),0,0,i).*...
                                GradJacobiP(XYZr(:,2),0,0,j).*...
                                    JacobiP(XYZr(:,3),0,0,k);
        GradChi_v(:,IndChi,3) =     JacobiP(XYZr(:,1),0,0,i).*...
                                    JacobiP(XYZr(:,2),0,0,j).*...
                                GradJacobiP(XYZr(:,3),0,0,k);
    end
    end
    end
end

return
