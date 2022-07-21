%%
% More documentation here:
%   help pLaplaceSolverNew
%   help NesterovPathFollowingNew
p = 1.0;   % The p of the p-Laplacian
% This next bit is a quick-and-dirty way of generating a FEM mesh
G = numgrid('S',52);
[I,J] = find(G);
V = [I J];
V = V-repmat(min(V),size(V,1),1);
V = V./repmat(max(V),size(V,1),1);
T = delaunay(V);
% Make up some boundary conditions
d = find(((V(:,2)>=0.25) & (V(:,2)<=0.75) & (V(:,1)==0)) | ((V(:,1)>0.6) & (V(:,2)>=0.25)));
g = zeros(size(V,1),1);
g(d) = 1;
% Zero forcing
f = zeros(size(T,1),1);
% Solve
[ u,SOL ] = pLaplaceSolverNew( T,V,f,p,g );
% Plot solution
close all;
disp(sprintf('done, its = %d, elapsed=%.3fs',length(SOL.lambda),SOL.elapsed));
patch('Faces',T,'Vertices',V,'FaceColor','interp','FaceVertexCData',u,'EdgeColor','none');
hold on;
plot(V(:,1),V(:,2),'.k','markers',0.5);
hold off;
colorbar;
axis square;
drawnow;
