function [Dx,Dy,W,M]=derivativesMatrix(f,v)

p=v(f(:,2),:)-v(f(:,1),:);
q=v(f(:,3),:)-v(f(:,1),:);
BT=[p q];
detB=(p(:,1).*q(:,2)-p(:,2).*q(:,1));
BTinv=[BT(:,4) -BT(:,2) -BT(:,3) BT(:,1)]./repmat(detB,1,4);
ni=[-1 -1;1 0;0 1];
m=size(p,1);
n=size(v,1);
Dx = sparse(m,n);
Dy = sparse(m,n);
BinvTu=@(u) [BTinv(:,1).*u(:,1)+BTinv(:,2).*u(:,2), ...
              BTinv(:,3).*u(:,1)+BTinv(:,4).*u(:,2)];
for i=1:3
    Grads = BinvTu(repmat(ni(i,:),m,1));
    Dx = Dx + sparse(1:m,f(:,i),Grads(:,1),m,n);
    Dy = Dy + sparse(1:m,f(:,i),Grads(:,2),m,n);
end
W = 0.5.*detB;
M = sparse(n,n);
MM = (1/24)*[2 1 1
             1 2 1
             1 1 2];
for i=1:3
    for j=1:3
        M = M + sparse(f(:,i),f(:,j),detB*MM(i,j),n,n);
    end
end
%M = spdiags(sum(M)',0,size(M,1),size(M,2));

