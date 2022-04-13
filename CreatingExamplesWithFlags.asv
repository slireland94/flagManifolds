%% 
% Seth Ireland, Anthony Caine and Brian Collery 
% 4.13.2022

% The purpose of this code is to create flags with different
% representatives to see where the shortest path occurs by testing them
% all. In this particular code we will test on full flags but the next
% example will be on flags that are not full flags but have a circular
% option to chose from. We will use a sampling method to try figure out
% which has the shortest distance. The hope is that this will lead to ideas
% about how to fix the code we have.

p =[1,1,1];
[l1,l2,FRepresentatives,flagLog,UFlag,SFlag,VFlag] = generateExamples(p);
function [l1,l2,FRepresentatives,flagLog,UFlag,SFlag,VFlag] =generateExamples(p)
    % For the moment we will assume one of our flags is the identity as it
    % is hopefully easier to see what is going on. So G = eye(n) which we
    % do not need.
    n = sum(p);
    F = specialOrtho(n);
    

    % Next we want to take different representatives of f

    [FRepresentatives] = generateExcessFi(F,n);

    % Now we want to see the shortest path by testing each representative.
    
    flagLength = lengthOfFlagRepresentatives(FRepresentatives,p,n);

    % See what happens if we use the log versions.
    flagLog = takeLogOfRepresentatives(FRepresentatives);
    flagLength2= lengthOfFlagRepresentatives(flagLog,p,n);
    % Next start plotting (mwahahahah) the graph of the distances.



    % plot( flaglength, the length of the flag)

    % Info to examine to see what is happening.
    l1 = flagLength;
    l2 = flagLength2;
    [UFlag,SFlag,VFlag] = svdFRep(FRepresentatives);
end



function [F] = specialOrtho(n)
A = rand(n);
R = zeros(n);
% implement Gram-Schmidt process to get a 'random' element of O(n)
for j = 1:n
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = F(:,i)' * A(:,j);
        v = v - R(i,j) * F(:,i);
    end
    R(j,j) = norm(v);
    F(:,j) = v / R(j,j);
end
F(:,1) = F(:,1)*det(F); % take Q\in O(n) and force it to be Q\in SO(n)
end


% get representatives of [Q] using cover of flag manifold
function [FiShort] = generateFi(F,p)
    d = length(p);
    pAlt = altSyntax(p);
    count_rows = height(F);
    count_columns = width(F);
    FiShort = zeros(count_rows,count_columns,2^(d-1));
    i = 1;
    FiShort(:,:,i) = F;
    nchoosek((0:3),3);
    columns = (1:count_columns);

    for j = 1:floor(d/2)
        C = nchoosek(pAlt,2*j);
            for k =1:size(C)
                i = i+1;
                FiShort(:,:,i) = F;
                FiShort(:,C(k,:),i) = -FiShort(:,C(k,:),i);
            end
    end
end

function [Fi] = generateExcessFi(F,n)
    % This function generates all possible Qi's rather than good
    % representatives.

    count_rows = height(F);
    count_columns = width(F);
    Fi   = zeros(count_rows,count_columns,2^(n-1));
    i = 1;
    Fi(:,:,i) = F;
    nchoosek((0:3),3);
    columns = (1:count_columns);

    for j = 1:floor(count_columns/2)
        C = nchoosek(columns,2*j);
            for k =1:size(C)
                i = i+1;
                Fi(:,:,i) = F;
                Fi(:,C(k,:),i) = -Fi(:,C(k,:),i);
            end
    end
end

function [pAlt] =altSyntax(p)
    pAlt = p;
    for i = 2:length(p)
        pAlt(i) = pAlt(i) + pAlt(i-1);
    end
end


% This function gives us the length of getting from the identity 
function [flagRepLength] = lengthOfFlagRepresentatives(FRepresentatives,p,n)
    flagRepLength = zeros(length(FRepresentatives),1);
    for i=1:length(FRepresentatives)
        flagRepLength(i) = flagMetric(FRepresentatives(:,:,i));
    end
end

function [l] = flagMetric(F)
    l = .5*trace(F'*F);

end

function [flagLog] = takeLogOfRepresentatives(FRepresentatives)
    k = size(FRepresentatives);
    flagLog = zeros(k(1),k(2),length(FRepresentatives));
    for i=1:length(FRepresentatives)
        flagLog(:,:,i) = logm(FRepresentatives(:,:,i));
    end
end

function [UFlag,SFlag,VFlag] = svdFRep(FRepresentatives)
    k = size(FRepresentatives);
    UFlag = zeros(k(1),k(2),length(FRepresentatives));
    SFlag = zeros(k(1),k(2),length(FRepresentatives));
    VFlag = zeros(k(1),k(2),length(FRepresentatives));
    for i=1:length(FRepresentatives)
        [U,S,V] = svd(FRepresentatives(:,:,i));
        UFlag(:,:,i) = U;
        SFlag(:,:,i) = S;
        VFlag(:,:,i) = V;
    end
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
