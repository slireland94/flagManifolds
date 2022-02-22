%% 
% Seth Ireland, Anthony Caine and Brian Collery 
% 1.26.2022

close all
clear all
clc

% This script is to create a geodesic on a Grassmannian with the potential
% to try do it on a flag in future. 

% what kind of Grassmannian do we have?
p = [2,2,6];
testCases = 4;
% Forgot the bar part
G1Bar = eye(10);
G2Bar = eye(10);

check =specialOrtho(testCases);
G2Bar(2:(1+testCases),2:(1+testCases)) = check;
G1Bar(2:(1+testCases),2:(1+testCases)) = specialOrtho(testCases);
[G1Bar,G2Bar] =barMatrix(G1Bar,G2Bar,p)
tic;
[G,A] = AltTestRun(p,G1Bar,G2Bar);

toc;
% testRun(p)
% in the Grassmannian case a_1 = sin^{-1}b_1


function [G,A,G1Bar,G2Bar] = AltTestRun(p,G1Bar,G2Bar)
    n = sum(p);
    % get our points in the Grassmannian (just representatives right now)
    % Thicken the bar of each grassmannian
    % respectively. I.e multiply them by their
    % SVD vector decomp respectively and then find their complement.
    G = createExpM(G1Bar,G2Bar,p,n);
    if(det(G) < 0.5)
        G(:,1) = -G(:,1);
        A = logm(G);
    else
        A = logm(G);
    end
end

function [G] =createExpM(G1Bar,G2Bar,p,n)
l = length(p);
G = zeros(n);
for i=1:n
    G(i,i) = G1Bar(:,i)'*G2Bar(:,i);
end
if( l > 2) %If SVD breaks down.
    pAlt = [0 altSyntax(p)];
    % Where it breaks down.
    for i=1:(l-1)
        % choice of A's
        for j = (pAlt(i)+1):pAlt(i+1)
            for k = (pAlt(i+1)+1):pAlt(l+1)
               G(j,k) = G1Bar(:,j)'*G2Bar(:,k);
               G(k,j) = G1Bar(:,k)'*G2Bar(:,j);
            end
        end
    end
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

    G1Slice = G1Bar(:,1:p(1));
    G2Slice = G2Bar(:,1:p(1));
    
    [U,S,V ] = svd(G1Slice'*G2Slice);
    G1Bar(:,1:p(1)) = G1Slice*U;
    G2Bar(:,1:p(1)) = G2Slice*V;

    for i =1:(l-1)
        G1Slice = G1Bar(:,(pAlt(i)+1):(pAlt(i+1)));
        G2Slice = G2Bar(:,(pAlt(i)+1):(pAlt(i+1)));
        [Uc,S,Vc ] = svd(G1Slice'*G2Slice);
        G1Bar(:,(pAlt(i)+1):(pAlt(i+1))) = G1Slice*Uc;
        G2Bar(:,(pAlt(i)+1):(pAlt(i+1))) = G2Slice*Vc;
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
function [Q] = specialOrtho(n)
A = rand(n);
R = zeros(n);
% implement Gram-Schmidt process to get a 'random' element of O(n)
for j = 1:n
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,i)' * A(:,j);
        v = v - R(i,j) * Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v / R(j,j);
end
Q(:,1) = Q(:,1)*det(Q); % take Q\in O(n) and force it to be Q\in SO(n)
end
