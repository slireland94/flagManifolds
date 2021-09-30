% Seth Ireland and Brian Collery 
% 9.7.2021

close all
clear all
clc

% this script replicates the iterative algorithm described in Ma's paper
% as of 9/13/2021, it does this for SPECIFIC REPRESENTATIVES Q1 and Q2
% a more general approach allowing for all representative pairs (oriented
% flag manifold) will come later

% what kind of flag do we have?
p = [2,2,2];
testRun(p);




function [dis] = testRun(p)
    n = sum(p);
    % get our points in the flag (just representatives right now)
    Q1 = specialOrtho(n);
    Q2 = specialOrtho(n);
    Q = Q1'*Q2;
    [H,G] = computeHG(Q,p);
    % Using the structure of the flag manifold
    d = length(p);
    [Qi] = generateQi(Q,p);
    % desiShort = zeros(2^(d-1),1);
    dis = zeros(2^(d-1),1);
    eigenValuesQ = zeros(2^(d-1),n);
    % Looking at the distance with the flag manifold structure
    for i= 1: 2^(d-1)
        eigenValuesQ(i,:) = eig(Qi(:,:,i));
        eigenValuesQ(i,:);
        [H,G] = computeHG(Qi(:,:,i),p);
        dis(i) =  sqrt(0.5*(trace(H'*H)));
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

function [Qi] = generateExcessQi(Q,n)
    % This function generates all possible Qi's rather than good
    % representatives.

    count_rows = height(Q);
    count_columns = width(Q);
    Qi   = zeros(count_rows,count_columns,2^(n-1));
    i = 1;
    Qi(:,:,i) = Q;
    nchoosek((0:3),3);
    columns = (1:count_columns);

    for j = 1:floor(count_columns/2)
        C = nchoosek(columns,2*j);
            for k =1:size(C)
                i = i+1;
                Qi(:,:,i) = Q;
                Qi(:,C(k,:),i) = -Qi(:,C(k,:),i);
            end
    end
end

function [pAlt] =altSyntax(p)
    pAlt = p;
    for i = 2:length(p)
        pAlt(i) = pAlt(i) + pAlt(i-1);
    end
end


function [H,G] = computeHG(Q,p)
    % initialize
    G0 = blockDiagSkewSym(p);
    H_hat = logm(Q*expm(G0)');
    H = projectToComp(H_hat,p);
    % run the algorithm
    error = 1;
    tolerance = 0.01;
    while error > tolerance  % tolerance from last to current
        if countNegEig(expm(H)'*Q,p) > 0 % if expm(H)'*Q has negative eigenvalues, then let's see it
            expm(H)'*Q
        end
        G_hat = logm(expm(H)'*Q);
        G = projectToWP(G_hat,p);
        if countNegEig(Q*expm(G)',p) > 0 % if Q*expm(G)' has negative eigenvalues, then let's see it
            Q*expm(G)'
        end
        H_hat = logm(Q*expm(G)');
        H = projectToComp(H_hat,p);
        errorM = Q - expm(H)*expm(G); %error matrix
        % 
        error = sqrt(0.5*trace(errorM'*errorM));
    end
    % distance = sqrt(0.5*trace(H'*H)) % distance between Q1 and Q2 ?
    % distance = sqrt(0.5*trace(G'*G)) % ???

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
        A = skewSym(p(j));
        B = blkdiag(B,A);
    end
end


function [G0] = skewSym(n)
R = rand(n);
G0 = R - R';
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
