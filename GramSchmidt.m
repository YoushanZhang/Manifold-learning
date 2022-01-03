function [R,Q]=GramSchmidt(A, Q)

 [m, n]=size(A);
 R=zeros(m,n);


 for j=1:1:n  
      v = A(:,j); 
   for i=1:1:j-1
      R(i, j) = dot(Q(:,i), A(:,j));
      v = v - R(i, j) * Q(:,i);
    end
      R(j, j) = norm(v, 'fro');
      Q(:,j) = v / R(j, j);
 end

end