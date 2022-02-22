%%Playing with geodesic on grassman
%Generating 2 points on Gr(2,4)
A=rand(4,2);
[q,r]=qr(A);
A=q(:,1:2); %finding an ON basis for col(A)

B=rand(4,2);
[q,r]=qr(B);
B=q(:,1:2); %finding an ON basis for col(A)

%% Finding ideal basis and principle angles
[u,s,v]=svd(A'*B);
A=A*u;
B=B*v;  %idea: A'B=usv' -> u'A'Bv=s
% https://math.stackexchange.com/questions/4167802/get-a-rotation-matrix-which-rotates-4d-vector-to-another-4d-vector?rq=1
%% building the geodesic
w=B-A*A'* B; % B minus its projection onto A.
w(:,1)=w(:,1)/norm(w(:,1)); %normalize difference vector
w(:,2)=w(:,2)/norm(w(:,2)); %normalize difference vector

dt=.01;
t=0:dt:1;

theta=acos(diag(s)); %principle angles as column vector

gamma=zeros([size(A) length(t)]); %allocating geodesic memory
for i=1:length(t)
    gamma(:,:,i)=A*diag(cos(t(i)*theta))+w*diag(sin(t(i)*theta));
end

