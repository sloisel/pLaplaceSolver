function [ u,SOL ] = pLaplaceSolverNew( T,V,fmid,p,g,varargin )
%PLAPLACESOLVERNEW Solves a p-Laplace type equation\
%   Solves a p-Laplace type equation L u = f.
%
%   Parameters:
%   T    The array of triangles (size is m by 3)
%   V    The array of vertices (size is n by 2)
%   fmid The forcing (this is a function handle)
%   p    The parameter of the p-Laplacian
%   g    The Dirichlet boundary conditions (size is n by 1)
%   Output:
%   u    The solution of the p-Laplace type equation.
%
%   Note: the mass matrix M is automatically computed from
%         T and V.
%
%   Optional parameters:
%   Any optional parameters are passed straight through to the
%   function NesterovPathFollowingNew. See
%
%   help NesterovPathFollowingNew

if(size(V,2)==2)
    [Dx,Dy,omega] = derivativesMatrix(T,V);
    Dz = sparse(size(Dx,1),size(Dx,2));
    D = {Dx,Dy,Dz};
else
    [D,omega]=stiffnessMatrix3d(T,V,1);
    Dx = D{1};
    Dy = D{2};
    Dz = D{3};
end
b = find_boundary(T);
i = setdiff(1:size(V,1),b);
W = spdiags(omega,0,length(omega),length(omega));
A = Dx'*W*Dx + Dy'*W*Dy;
AII = A(i,i);
AIB = A(i,b);
g(i) = AII\(-AIB*(g(b)));
ni = length(i);
m = size(Dx,1);
Dxi = Dx(:,i);
Dyi = Dy(:,i);
Dzi = Dz(:,i);
if(size(V,2)==2)
    Di = {Dxi,Dyi};
else
    Di = {Dxi,Dyi,Dzi};
end
dgdx = Dx*g;
dgdy = Dy*g;
dgdz = Dz*g;
pinf = false;
preal = p;
if(p==inf)
    pinf = true;
    q = 1;
    gnormp = max(dgdx.^2+dgdy.^2+dgdz.^2).^(1/2);
    Rstar = 2*(1+gnormp);
    fnorm = omega'*abs(fmid);
    p = 1;
else
    q = p/(p-1);
    if(q==inf)
        fnorm = max(abs(fmid));
    else
        fnorm = (omega'*abs(fmid).^q).^(1/q);
    end
    gnormp = omega'*((dgdx.^2+dgdy.^2+dgdz.^2).^(p/2));
    Rstar = 2*(1+gnormp);
end
Mmid = sparse(repmat(1:m,1,size(T,2))',T(:),1/size(T,2),m,size(V,1));
MmidI = Mmid(:,i);

L = max(V(:,1))-min(V(:,1));
if(p==1)
    assert(L*fnorm<1,sprintf('L*fnorm=%f>=1',L*fnorm));
    R = (2+2*gnormp/(1-L*fnorm));
    if(pinf), R = R*max(omega); end
else
    R = max(2 + 4*(p/2)^(1/(1-p))*(p-1)*L^q*fnorm^q+8*gnormp);
end
if(pinf)
    s = max(1+((dgdx).^2+(dgdy).^2+(dgdz).^2).^(1/2));
    C = 1;
else
    s = 1+((dgdx).^2+(dgdy).^2+(dgdz).^2).^(p/2);
    C = omega;
end
assert(max(omega.*s)<R);
F = pLaplaceBarrierY( preal,Di,D,g,R,omega);
x = [zeros(ni,1);s];
ff = -Mmid'*(omega.*fmid);
ffI = -MmidI'*(omega.*fmid);
[ x,SOL ] = NesterovPathFollowingNew( [ffI;C],F,x,varargin{:} );
u = g;
u(i) = u(i) + x(1:ni);
if(pinf)
    SOL.J = @(u) max((((Dx*u+dgdx).^2+(Dy*u+dgdy).^2+(Dz*u+dgdx).^2).^0.5))+ff'*u;
    SOL.unormp = max(omega)*max((((Dx*u).^2+(Dy*u).^2+(Dz*u).^2).^0.5));
else
    SOL.J = @(u) (1/p)*(omega'*(((Dx*u+dgdx).^2+(Dy*u+dgdy).^2+(Dz*u+dgdx).^2).^(p/2)))+ff'*u;
    SOL.unormp = omega'*(((Dx*u).^2+(Dy*u).^2+(Dz*u).^2).^(p/2));
end
SOL.R = R;
SOL.u = u;
SOL.fmid = fmid;
SOL.g = g;
SOL.T = T;
SOL.V = V;
SOL.p = p;
SOL.gnormp = gnormp;
assert(SOL.unormp<SOL.R/2,'Something horrible happened and SOL.unormp is larger than SOL.R/2');
end
