function indLogic=upperLogic(p)
indmx = reshape([1:p^2],p,p); 
ind = indmx(triu(indmx,1)>0);  % upperInd
indLogic=zeros(p);
indLogic(ind)=1;
indLogic=logical(indLogic);