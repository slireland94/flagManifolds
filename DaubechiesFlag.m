%% 
% Seth Ireland, Anthony Caine and Brian Collery 



close all
clear all
clc


% if k>n, subspace becomes dependent on circshift function, Matlab doesn't compute daub(k,1) for k>49
n = 14; % computations will take place on Gr(n,2n)
k = [1,5]; % the 2 k values we want to measure distance between
p = [n,n];
t = 0:0.01:1;
% makeD2k generates a representative for D_{2k} as a point on Gr(n,2n)
%makeD2k(k,n)
G1 = makeD2k(k(1),n);
G2 = makeD2k(k(2),n);
% G1 = flagDaubechies(G1,p);
% G2 = flagDaubechies(G2,p);
[A,G1Bar,G2Bar] = geodesicMat(G1,G2,p);
I= imread('MNISTExample.png');

Im = im2double(I);
ImPlots = zeros(28,28,length(t));
ImPlotsFlag = zeros(28,28,length(t));
ImPlotsRand = zeros(28,28,length(t));
B = specialOrtho(28);
Blog = logm(B);
for i=1:9
    
    expA =(G1Bar * expm(A*(t(i))));
    expA = expA(:,(1:14));
    expB = (G1Bar * expm(Blog*(t(i))));
    expB = expB(:,(1:14));
    projB = expB*expB';
    projA = expA*expA';
    expFlag = expA(:,1:7);
    projFlag = expFlag*expFlag';
    ImPlots(:,:,i) = projA*Im*projA;
    ImPlotsFlag(:,:,i) = projFlag*Im*projFlag;
    ImPlotsRand(:,:,i) = projB*Im*projB;
    frameName = strcat('frame00',int2str(i),'.png');
    saveas(imshow(imresize(ImPlotsFlag(:,:,i),13)),frameName);
end
for i=10:(length(t)-2)
    
    expA =(G1Bar * expm(A*(t(i))));
    expA = expA(:,(1:14));
    expB = (G1Bar * expm(Blog*(t(i))));
    expB = expB(:,(1:14));
    projB = expB*expB';
    projA = expA*expA';
    expFlag = expA(:,1:7);
    projFlag = expFlag*expFlag';
    ImPlots(:,:,i) = projA*Im*projA;
    ImPlotsFlag(:,:,i) = projFlag*Im*projFlag;
    ImPlotsRand(:,:,i) = projB*Im*projB;
    frameName = strcat('frame0',int2str(i),'.png')
    saveas(imshow(imresize(ImPlotsFlag(:,:,i),13)),frameName);
end
for i=(length(t)-1):(length(t))
    
    expA =(G1Bar * expm(A*(t(i))));
    expA = expA(:,(1:14));
    expB = (G1Bar * expm(Blog*(t(i))));
    expB = expB(:,(1:14));
    projB = expB*expB';
    projA = expA*expA';
    expFlag = expA(:,1:7);
    projFlag = expFlag*expFlag';
    ImPlots(:,:,i) = projA*Im*projA;
    ImPlotsFlag(:,:,i) = projFlag*Im*projFlag;
    ImPlotsRand(:,:,i) = projB*Im*projB;
    frameName = strcat('frame',int2str(i),'.png')
    saveas(imshow(imresize(ImPlotsFlag(:,:,i),13)),frameName);
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
    % Thicken the bar of each grassmannian
    % respectively. I.e multiply them by their
    % SVD vector decomp respectively and then find their complement.
    G1Bar = Thicken(G1,p);
    G2Bar = Thicken(G2,p);
    [G1Bar,G2Bar] = barMatrix(G1Bar,G2Bar,p);
    Q = G1Bar'*G2Bar;
    [A,G] = computeHG(Q,p);
end

function [G] = flagDaubechies(G,p)
    tempG = G;
    sep = p(1);
    for i=1:p(1)
        G(:,i) = (1/sqrt(2))*(tempG(:,2*(i-1) + 1) + tempG(:,2*(i-1) + 2));
    end
    for i=1:p(1)
        G(:,i+sep) = (1/sqrt(2))*(tempG(:,2*(i-1) + 1) - tempG(:,2*(i-1) + 2));
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



function [H,G,count] = computeHG(Q,p)
    % initialize
    origQ = Q;
    countMax = 50000;
    G0 = blockDiagSkewSym(p);
    H_hat = logOfMatrix(Q*expm(G0)');
    H = projectToComp(H_hat,p);
    % run the algorithm
    error = 1;
    tolerance = 0.01;
    count=[0,0];
    while (error > tolerance) && (count(1) <countMax)  % tolerance from last to current
        count(1)=count(1)+1;
        if countNegEig(expm(H)'*Q,p) > 0 % if expm(H)'*Q has negative eigenvalues, then let's see it
            writematrix(expm(H)'*Q,'examples.xls');
            writematrix(p,'pUsed.xls');
            disp('got here');
            writematrix([1],'decision.xls');
            count(2)=count(2)+1
        end
        G_hat = logOfMatrix(expm(H)'*Q);
        if ( (max(max(abs(expm(G_hat)-expm(H)'*Q)))) > .000001 )
            G_hat
            expm(H)'*Q
            disp('Got Here')
            disp(max(max(abs(expm(G_hat)-expm(H)'*Q))))
            pause
        end
        G = projectToWP(G_hat,p);
        if countNegEig(Q*expm(G)',p) > 0 % if Q*expm(G)' has negative eigenvalues, then let's see it
            writematrix(Q,'examples.xls');
            writematrix(p,'pUsed.xls');
            writematrix([2],'decision.xls');
        end
        H_hat = logOfMatrix(Q*expm(G)');
        if ( (max(max(abs(expm(H_hat)-Q*expm(G)')))) > .000001 )
            H_hat
            disp('got here')
            disp(max(max(abs(expm(H_hat)-Q*expm(G)'))))
            Q*expm(G)'
            pause
        end
        H = projectToComp(H_hat,p);
        % Make the error less than max of max
        
        errorM = Q - expm(H)*expm(G); %error matrix 
        error = max(max(abs(errorM)));
    end
    if (count(1) == 1000)
        disp('Found a Q')
        writematrix(origQ,'examples.xls');
        writematrix(p,'pUsed.xls');
    end
    if (count(1) == countMax)
        disp('Reached max iterations.')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G] = projectToWP (G_hat,p)
    A = zeros(sum(p));
    topLeftCorner=1;
    for j = 1:length(p)
        blockSize = p(j);
        for k = 0:blockSize-1
            for l = 0:blockSize-1
                A(topLeftCorner+k,topLeftCorner+l) = G_hat(topLeftCorner+k,topLeftCorner+l);
            end
        end 
        topLeftCorner = p(j) + topLeftCorner;
    end
    G = A;
end

function [H] = projectToComp (H_hat,p)
    topLeftCorner=1;
    for j = 1:length(p)
        blockSize = p(j);
        for k = 0:blockSize-1
            for l = 0:blockSize-1
                H_hat(topLeftCorner+k,topLeftCorner+l)=0;
            end
        end 
        topLeftCorner = p(j) + topLeftCorner;
    end
    H = H_hat;
end

function [B] = blockDiagSkewSym(p)
    B =[];
    for j = 1:length(p)
        A = skewSymRand(p(j));
        B = blkdiag(B,A);
    end
end

function [mat] = logOfMatrix(M)
    mat = logm(M);
    mat = skewMatrix(mat);
    
end

function [G0] = skewSymRand(n)
R = rand(n);
G0 = R - R';
end

function [skewM] =skewMatrix(M)
    skewM = .5*(M - M');
end



% check if any eigenvalues of a matrix are negative
function [TF] = countNegEig(B,p)
    l = length(B);
    E = eig(B);
    TF = 0;
    for j = 1:l
        if imag(E(j)) == 0 && real(E(j)) < 0
            TF = TF + 1; % B has a negative eigenvalue
        else
            TF = TF; % B has no negative eigenvalues
        end
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

