function [T,V,R] = subdivide(T,V)
% [T,V,R] = subdivide(T,V) -- Subdivide a mesh
%
% Example:
% [T,V] = subdivide([1 2 3;1 3 4],[0 0;1 0;1 1;0 1])
% x = V(:,1);
% y = V(:,2);
% patch(x(T'),y(T'),rand(size(T')));
x = V(:,1);
y = V(:,2);
n = length(x);
edges = [T(:,1) T(:,2)
         T(:,2) T(:,3)
         T(:,1) T(:,3)];
edges = [min(edges,[],2) max(edges,[],2)];
m = size(T,1);
m1 = (n+1:n+m)';
m2 = (n+m+1:n+2*m)';
m3 = (n+2*m+1:n+3*m)';
T2 = [ T(:,1) m1 m3
       T(:,2) m2 m1
       T(:,3) m3 m2
       m1     m2 m3 ];
midx = 0.5*(x(edges(:,1)) + x(edges(:,2)));
midy = 0.5*(y(edges(:,1)) + y(edges(:,2)));
x2 = [x;midx];
y2 = [y;midy];
[b,r,s] = unique(edges,'rows');
r = [(1:n)' ; r+n];
s = [(1:n)' ; s+n];
x2 = x2(r);
y2 = y2(r);
T2 = s(T2);
N = size(T,1);
R = max(max(sparse(1:length(x),1:length(x),1,length(x),length(x2)),...
            sparse(T2(1:3*N,1),T2(1:3*N,2),0.5,length(x),length(x2))),...
        sparse(T2(1:3*N,1),T2(1:3*N,3),0.5,length(x),length(x2)));
T = T2;
V = [x2 y2];

