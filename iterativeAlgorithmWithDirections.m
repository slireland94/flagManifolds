% Anthony Caine, Seth Ireland and Brian Collery 
% 7.18..2022

close all
clear all
clc

% this script replicates the iterative algorithm described in Ma's paper
% as of 9/13/2021, it does this for SPECIFIC REPRESENTATIVES Q1 and Q2
% but it allows to have different partitions.

% what kind of flag do we have? The partition is what sets of vectors are
% grouped together to form the partially orientated flag manifold. I used
% cell notation as we can have elements of our partition of different
% length and that was causing issues with the matrix.
p = [3,3,3,3,3,3,3];
partition = {[1,3,5],[2,4],[6,7]};
tic
[dis,H,Qi,compare,amountOfReps]=testRun(p,partition);
toc



function [dis,tangents,Qi,compare,amountOfReps] = testRun(p,partition)
    n = sum(p);
    d = length(p);
    amountOfReps = partAmount(partition);
    % get our points in the flag (just representatives right now)
    Q1 = eye(n);
    Q2 = specialOrtho(n);
%     Q1 = specialOrtho(n)
%     Q2 = specialOrtho(n)
    Q = Q1'*Q2;
    [H,G] = computeHG(Q,p);
    tangents = zeros(n,n,amountOfReps);
    compare = zeros(n,n,amountOfReps);
    % Using the structure of the flag manifold
    d = length(p);
    [Qi] = generateQiDirection(Q,p,partition,amountOfReps);
    dis = zeros(amountOfReps,1);
    eigenValuesQ = zeros(amountOfReps,n);
    % Looking at the distance with the flag manifold structure
    for i= 1: amountOfReps
        [H,G] = computeHG(Qi(:,:,i),p);
        dis(i) =  sqrt(0.5*(trace(H'*H)));
        tangents(:,:,i)= H;
        compare(:,:,i) = Qi(:,:,i)'*expm(H);
    end
    eigenValuesQ;
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

    


% get representatives of [Q] using cover of flag manifold
function [QiShort] = generateQi(Q,p)
    d = length(p);
    pAlt = altSyntax(p);
    count_rows = height(Q);
    count_columns = width(Q);
    QiShort = zeros(count_rows,count_columns,2^(d-1));
    i = 1;
    QiShort(:,:,i) = Q;
    nchoosek((0:3),3);
    columns = (1:count_columns);

    for j = 1:floor(d/2)
        C = nchoosek(pAlt,2*j);
            for k =1:size(C)
                i = i+1;
                QiShort(:,:,i) = Q;
                QiShort(:,C(k,:),i) = -QiShort(:,C(k,:),i);
            end
    end
end

function [QiShort] = generateQiDirection(Q,p,partition,amountOfReps)
    d = length(p);
    pAlt = altSyntax(p);
    count_rows = height(Q);
    count_columns = width(Q);
    QiShort = zeros(count_rows,count_columns,amountOfReps);
    i = 1;
    l = 1;
    QiShort(:,:,i) = Q;
    partLength = length(partition);
    % Here we are going to have to use the first partition to create the
    % first set of cases. The reason the first partition is created in a
    % different way is that it can have the identity and the rest cannot.
    
    
    % We need a for loop to go through each partition
    for m= 1:partLength
        % This is the part where we stack representations
        d =length(partition{m})
        for j = 1:floor(d/2)
            pChoice = zeros(length( partition{m}),1);
            for q=1:length( partition{m})
                pChoice(q) =pAlt(partition{m}(q) );
            end
            C = nchoosek( pChoice,2*j);
            for s = 1:i
                for k =1:size(C)
                    l = l+1;
                    QiShort(:,:,l) = QiShort(:,:,s);
                    QiShort(:,C(k,:),l) = -QiShort(:,C(k,:),l);
                end
            end
        end
        % This is the flattening part
        i=l
    end
    l
end

function [pAlt] =altSyntax(p)
    pAlt = p;
    for i = 2:length(p)
        pAlt(i) = pAlt(i) + pAlt(i-1);
    end
end


function [H,G] = computeHG(Q,p)
    % initialize
    origQ = Q;
    G0 = blockDiagSkewSym(p);
    H_hat = logOfMatrix(Q*expm(G0)');
    H = projectToComp(H_hat,p);
    % run the algorithm
    error = 1;
    tolerance = 0.0001;
    count=[0,0];
    while (error > tolerance) && (count(1) <10000)  % tolerance from last to current
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
    
end

% Questions:
% how to compute 'unique singular values' of H
% distance = sqrt(lambda_1^2 + lambda_2^2)

% checkLog checks to see if Q*expm(G0)' or expm(H)'*Q ever have negative
% eigenvalues


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


function [partitionPossibilities] = partAmount(partition)
    partitionPossibilities = 1;
    for i = 1:length(partition)
        partitionPossibilities = partitionPossibilities*2^(length(partition{i})-1);
    end
end
