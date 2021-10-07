% This file is to examine matrices where the log function on its projection
% will not work. One issue is that we do not have p from the text file so
% we will have to alter that manually. 
p = [2,2,2];
A = readmatrix('examples.xls');
eig(A)
[eigenValues] = piecewiseProjection(p,A)
% First way to examine it is to check how its eigenvalues change when we do
% "piecewise" projection. This is when we do the projection for each step
% but check its eigenvalues in between and see if we can find a pattern.



function [eigenValues] = piecewiseProjection(p,Q)
    d = length(p);
    n = sum(p);
    eigenValues = zeros(n,(d+1));
    eigenValues(:,1) = eig(Q);
    topLeftCorner=1;
    for j = 1:d
        blockSize = p(j);
        for k = 0:blockSize-1
            for l = 0:blockSize-1
                Q(topLeftCorner+k,topLeftCorner+l)=0;
            end
        end 
        eigenValues(:,j+1) = eig(Q);
        topLeftCorner = p(j) + topLeftCorner;
        
    end

end