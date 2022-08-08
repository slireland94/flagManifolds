%% 
% Seth Ireland, Anthony Caine and Brian Collery 



close all
clear all
clc


% if k>n, subspace becomes dependent on circshift function, Matlab doesn't compute daub(k,1) for k>49
n = 4; % computations will take place on Gr(n,2n)
k = [1,2]; % the 2 k values we want to measure distance between
p = [n,n];
t = 0:0.01:5;
% makeD2k generates a representative for D_{2k} as a point on Gr(n,2n)
%makeD2k(k,n)
G1 = makeD2k(k(1),n);

G2 = makeD2k(k(2),n);
[A,G1Bar,G2Bar] = geodesicMat(G1,G2,p);
I= imread('MNISTExample.png');

Im = im2double(I);
ImPlots = zeros(2*n,2*n,length(t));
B = rand(28,14);
for i=1:length(t)
    
    expA =(G1Bar * expm(A*(t(i))));
    expA = expA(:,(1:n));
    projA = expA*expA';
%     ImPlots(:,:,i) = projA*Im*projA;
%     imshow(imresize(ImPlots(:,:,i),1))
end



function [Q] = makeD2k(k,n)
    D = dauboMat(k,n);
    Q = gramSchmidt(D,n);
end


function [D] = dauboMat(k,n)
    d = dbaux(k,1);
    v = zeros(2*n,1);
    v(1:2*k)=d'; % this is the basic column...now need to circshift it
    for i=1:n
        D(:,i) = circshift(v,2*i-2); % no matter where we start the circshift, we get same subspace (as long as k<=n)
    end
    D = D(1:2*n,:);
end

function [Q] = gramSchmidt(A,n)
for j = 1:n
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,i)' * A(:,j);
        v = v - R(i,j) * Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v / R(j,j);
end
%Q(:,n)=det(Q)*Q(:,n);
end


function [A,G1Bar,G2Bar] = geodesicMat(G1,G2,p)
    n = sum(p);
    % Get the SVD decom of their product.
    B = G1'*G2;
    [U,S,V] = svd(B);
    % Thicken the bar of each grassmannian
    % respectively. I.e multiply them by their
    % SVD vector decomp respectively and then find their complement.
    G1Bar = Thicken(G1,p);
    G2Bar = Thicken(G2,p);
    [G1Bar,G2Bar] = barMatrix(G1Bar,G2Bar,p);
    G = G1Bar'*G2Bar;
    if (det(G) < 0)
        G2Bar(:,1) = G2Bar(:,1)*(-1);
        G = G1Bar'*G2Bar;
    end
    A = logm(G);
end



function [GThick] = Thicken(G,p)
    [WComp] = ComplementSpace(G, p);
    [GThick] = Fuse(G,WComp);
end


function [G1Bar,G2Bar] =barMatrix(G1Bar,G2Bar,p)
    l = length(p);
    pAlt =altSyntax(p);
    % We need to add the bar to the first section in a different manner
    % than the other sections as there is no zeroth element in a vector in
    % matlab. This is more of a coding issue than a math issue.

    % We will bar the first section by itself

    G1Slice = G1Bar(:,1:p(1));
    G2Slice = G2Bar(:,1:p(1));
    
    [U,S,V ] = svd(G1Slice'*G2Slice);
    G1Bar(:,1:p(1)) = G1Slice*U;
    G2Bar(:,1:p(1)) = G2Slice*V;

    for i =1:(l-1)
        G1Slice = G1Bar(:,(pAlt(i)+1):(pAlt(i+1)));
        G2Slice = G2Bar(:,(pAlt(i)+1):(pAlt(i+1)));
        [U,S,V ] = svd(G1Slice'*G2Slice);
        G1Bar(:,(pAlt(i)+1):(pAlt(i+1))) = G1Slice*U;
        G2Bar(:,(pAlt(i)+1):(pAlt(i+1))) = G2Slice*V;
    end
end


function [WComp] = ComplementSpace(W,p)
    % Appropriate Gram-Schmidt to this question.
    % Note that this is currently for a Grassmannian
    % but there is the bones to generalise to a flag manifold.
    l = length(p);
    n = sum(p);
    WComp = rand(n,p(l));
    
    % implement Gram-Schmidt process to orthogonalise the complement space
    % to the rest of 
    for j = 1:p(l)
        v = WComp(:,j);
        
        % This makes the new vector perpendicular to the old subspace.
        for i = 1:sum(p(1:(l-1)))
            const = W(:,i)' * WComp(:,j);
            v = v - const * W(:,i);
        end
        
        % This makes the new vector perpendicular to the new subspace.
        for i = 1:j-1
            const = WComp(:,i)' * WComp(:,j);
            v = v -const* WComp(:,i);
        end
        WComp(:,j) = v / norm(v);
    end

end

function [GBar] = Fuse(W,WComp,p)
% Here we fuse the two matrices together to make the Grassmannian. This
% will allow us to get our GBar fully. This is split as a function for when
% we get to flag manifolds.
GBar = [W WComp];
end

function [G] = Grassmannian(p)

% Note that this Grassmannian is not neccessarily in SO(n). It could be in
% O(n) but we need to view the complement space first i think?
l = length(p);
n = sum(p);
m = sum(p(1:(l-1)));
A = rand(n,m);
for j = 1:m
    v = A(:,j);
    for i = 1:j-1
        const = G(:,i)' * A(:,j);
        v = v - const * G(:,i);
    end
    R(j,j) = norm(v);
    G(:,j) = v / R(j,j);
end
%G(:,1) = G(:,1)*det(G); % take Q\in O(n) and force it to be Q\in SO(n)
end

function [pAlt] =altSyntax(p)
    pAlt = p;
    for i = 2:length(p)
        pAlt(i) = pAlt(i) + pAlt(i-1);
    end
end
