% getHDs.m - Obtaining Historical Decomposition
function HD = getHDs(IRF,n,u,whichVar)

eye_n = eye(n);
hmax = size(u,1)-1;
HD = nan(n,hmax+1);

for h = 0:hmax
   for j = 1:n
   HD(j,h+1) = eye_n(:,whichVar)'*IRF(:,:,h+1)*eye_n(:,j)*eye_n(:,j)'*u(end-h,:)'; 
   end    
end

HD = sum(HD,2);

end
