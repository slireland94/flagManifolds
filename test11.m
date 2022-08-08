% Seth Ireland, Anthony Caine and Brian Collery 
% 4.14.2022

close all
clear all
clc

% this script replicates the iterative algorithm described in Ma's paper
% as of 9/13/2021, it does this for SPECIFIC REPRESENTATIVES Q1 and Q2
% a more general approach allowing for all representative pairs (oriented
% flag manifold) will come later

% what kind of flag do we have?
p = [5,5,5];
n= sum(p);

% testRun uses Ma's algorithm on every possible representation of the flag.
tic
[dis,H,Qi,compare,Q1,Q2]=testRun(p);
toc


tic
% Here I was just trying to see the difference between the max and min
% representations of the transition map between the two flags.
 [M,I] = min(dis);
 [M,InMax]= max(dis);
% minMat = Qi(:,:,1)'*Qi(:,:,I);
% maxMat = Qi(:,:,1)'*Qi(:,:,InMax);
% [HMin,GMin] = computeHG(Qi(:,:,I),p);
% [HMax,GMax] = computeHG(Qi(:,:,InMax),p);


% This is to apply Ma's algorithm to a the representation that went through
% 'flag' principle angles.
[G1Bar,G2Bar,U,V] =barMatrix(eye(n),Qi(:,:,1),p);
G = G2Bar'*G1Bar;
tan = logm(G);
[H,Gsvd,count] = computeHG(G,p);


% This is just to track the different distances.
dmin = dis(I);
dmax = dis(InMax);
d2 = sqrt(.5*trace(H'*H));
diff = d2-dmin
toc

% tic
% [HU,GU] = computeHG(U,p);
% [HV,GV] = computeHG(V,p);
% toc

% Make the log shew symmetric real.


function [dis,tangents,Qi,compare,Q1,Q2] = testRun(p)
    n = sum(p);
    d = length(p);
    
    % get our points in the flag (just representatives right now)
    Q1 = eye(n);
    Q2 = specialOrtho(n);
%     Q1 = specialOrtho(n)
%     Q2 = specialOrtho(n)
    Q = Q1'*Q2;
    [H,G] = computeHG(Q,p);
    tangents = zeros(n,n,2^(d-1));
    compare = zeros(n,n,2^(d-1));
    % Using the structure of the flag manifold
    d = length(p);
    [Qi] = generateQi(Q,p);
    dis = zeros(2^(d-1),1);
    eigenValuesQ = zeros(2^(d-1),n);
    
    
    % Looking at the distance with the flag manifold structure
    for i= 1: 2^(d-1)
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


function [H,G,count] = computeHG(Q,p)
    % initialize
    origQ = Q;
    G0 = blockDiagSkewSym(p);
    H_hat = logOfMatrix(Q*expm(G0)');
    H = projectToComp(H_hat,p);
    % run the algorithm
    error = 1;
    tolerance = 0.00001;
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


function [G1Bar,G2Bar,U,V] =barMatrix(G1Bar,G2Bar,p)
    l = length(p);
    pAlt =altSyntax(p);
    % We need to add the bar to the first section in a different manner
    % than the other sections as there is no zeroth element in a vector in
    % matlab. This is more of a coding issue than a math issue.

    % We will bar the first section by itself
    G1 = G1Bar;
    G2 = G2Bar;
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
    U = G1'*G1Bar;
    V = G2'*G2Bar;
end



function [H,count] = computeUHV(Q,Q1,Q2,p)
    % initialize
    origQ = Q;
    [U,V] = barMatrix(Q1,Q2,p);
    tolerance = .001;
    count = 0;
    error = 100;
    while ((error > tolerance) && (count < 10000))
        H_Hat = logm(U'*Q*V);
        H = projectToComp(H_Hat,p);

        U_Hat = logm( Q*V*expm(H)');
        U_log = projectToWP(U_Hat,p);
        
        

        V_Hat = logm(Q'*U*expm(H));
        V_log = projectToWP(V_Hat,p);

        V = expm(V_log);

        U = expm(U_log);

        errorM = Q - U*expm(H)*V'; %error matrix 
        error = max(max(abs(errorM)));
        count = count +1;
    end
end
