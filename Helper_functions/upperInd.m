function ind=upperInd(p)
indmx = reshape([1:p^2],p,p); 
ind = indmx(triu(indmx,1)>0);