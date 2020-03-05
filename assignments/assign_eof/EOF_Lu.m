 function [PC,EOF,lmd,G,cumG,trans]=EOF_Lu(X,K)
% Solving X=EOF*PC.
%  
% X     : [m,n] space-by-time
% PC    : [K,n] principal components
% K     : number of modes you want to get
% lmd   : [K,1] eigenvalues in the descending order
% EOF   : [m,K] eigenvectors
% G     : [K,1] contribution ratio
% cumG  : [K,1] culmulative contribution
% trans : If space>>time, do spatio-temporal transformation!

%%   Tick out Nan     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete the bin that ever has NaN value (retrieved later). 
NanMatrics=isnan(X);       % isNaN matrix
NanVec=any(NanMatrics,2);  % Any NaN value along the time? space-by-1
NanIndx=find(NanVec==1);   % Spatial indices of bins that has NaN
% NoNIndx=find(NanVec==0);   % Spatial indices of bins that has no NaN
X(NanIndx,:)=[];           % Delete the 'nan' bins

%%   Covariance matrix      
[m,n]=size(X); % note the m here is not the real number of spatial grids
% If m>>n, do spatio-temporal transformation to save eig's time.
if m>10*n                      
    C=X'*X;                   % n-by-n
    trans=1;
else
    C=X*X';                   % m-by-m 
    trans=0;  
end
% Record C's size
mn=size(C,1);
%%   Eigenvalues and eigenvectors   
% C*V=V*D.  C,D~[x]^2, V~dimensionless (mn-by-mn)
[V,D]=eig(C,'balance');      
% clear C;
% put eigenvalues in descending order
[lmd,ind]=sort(diag(D),'descend'); 
% reorder columns of V
V=V(:,ind);             

%%   Contribution ratio            
G(1:K)=lmd(1:K)/sum(lmd);           % contribution from the first K modes
cumG=cumsum(G);  

%%  If transformed, modify the eigenvectors
% Note that eigenvalues of XX' and X'X are the same (first several modes), 
% while the eigenvectors not.
% Make V's size always m-by-mn 
if trans==1                           
    V=X*V/sqrt(diag(lmd));    % If m>>n, V is m-by-n; else V is m-by-m.
end                            

%%   Solving X = EOF * PC for PC
PC=V'*X;                      % PC~[x]:mn-by-n   

%%  Transfer the units
V=V*diag(sqrt(lmd));          % V~[x]
PC=diag(1./sqrt(lmd))*PC;     % make PC dimensionless 

%%   Retrieve NaN       
% m here is the total number of spatial grids
m=m+length(NanIndx);
EOF=zeros(m,mn);  
EOF(NanIndx,:)=NaN;
EOF(NanVec==0,:)=V;

%%  First K modes    
EOF=EOF(:,1:K);                           
lmd=lmd(1:K); 
PC=PC(1:K,:);
disp('%----------      EOF DONE    ------------%')

end 
