% This file is to examine matrices where the log function on its projection
% will not work. One issue is that we do not have p from the text file so
% we will have to alter that manually. 
clear
p = readmatrix('pUsed.xls');
A = readmatrix('examples.xls');
decide = readmatrix('decision.xls');
eig(A);
p
% Look over the code to see if it is actually orthogonal;
det(A)
% Look at the order off projection matrices.
A
computeHG(A,p)

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



function [H,G] = computeHG(Q,p)
    % initialize
    origQ = Q;
    G0 = blockDiagSkewSym(p);
    H_hat = logOfMatrix(Q*expm(G0)');
    H = projectToComp(H_hat,p);
    % run the algorithm
    error = 1;
    tolerance = 0.001;
    count=[0,0];
    while (error > tolerance) && (count(1) <200)  % tolerance from last to current
        count(1)=count(1)+1
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
        H
        G
        error
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

