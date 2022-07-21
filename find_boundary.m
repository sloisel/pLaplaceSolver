function [boundary,fboundary] = find_boundary3d( T )
%FUND_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
M = size(T,1);
N = size(T,2);
faces = zeros(M*N,N-1);
for k=1:size(T,2)
    faces((k-1)*M+1:k*M,:) = T(:,setdiff(1:N,k));
end
faces = sort(faces,2);
[b,m,n] = unique(faces,'rows');
counts = accumarray(n,1);
fboundary = b(counts==1,:);
boundary = unique(fboundary);
end

