%% 
% Seth Ireland, Anthony Caine and Brian Collery 
% 1.26.2022

close all
clear all
clc

% This script is to create a geodesic on a Grassmannian with the potential
% to try do it on a flag in future. 

% what kind of Grassmannian do we have?
p = [4,4];
[G,A,G1Bar,G2Bar]=testRun(p);

% testRun(p)


function [G,A,G1Bar,G2Bar] = testRun(p)
    n = sum(p);
    % get our points in the Grassmannian (just representatives right now)
    G1 = Grassmannian(p);
    G2 = Grassmannian(p);

    % Thicken the bar of each grassmannian
    % respectively. I.e multiply them by their
    % SVD vector decomp respectively.
    [G1Bar,G2Bar] = barMatrix(G1,G2,p);
    G = G1Bar'*G2Bar;
    if (det(G) < 0)
        G2Bar(:,1) = G2Bar(:,1)*(-1);
        G = G1Bar'*G2Bar;
    end
    D = det(G)
    A = logm(G);

    t = 0:0.01:1;
    for i=1:length(t)
        expA =expm(A*(t(i)));
        (G1Bar*expm(A*(t(i))))'*(G1Bar*expm(A*(t(i))));
    end
    
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
    G = G1Bar;
    G2 = G2Bar;
    G1Slice = G1Bar(:,1:p(1));
    G2Slice = G2Bar(:,1:p(1));
    
    [U,S,V ] = svd(G1Slice'*G2Slice);
    G1Bar(:,1:p(1)) = G1Slice*U;
    G2Bar(:,1:p(1)) = G2Slice*V;

    for i =1:(l-2)
        G1Slice = G1Bar(:,(pAlt(i)+1):(pAlt(i+1)));
        G2Slice = G2Bar(:,(pAlt(i)+1):(pAlt(i+1)));
        [Uc,S,Vc ] = svd(G1Slice'*G2Slice);
        G1Bar(:,(pAlt(i)+1):(pAlt(i+1))) = G1Slice*Uc;
        G2Bar(:,(pAlt(i)+1):(pAlt(i+1))) = G2Slice*Vc;
        i+1
    end
    U = G'*G1Bar
    V = G2'*G2Bar
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

