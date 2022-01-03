function [M, N]=randOrthoMat(A)

%   ret=0;
  while(det(A) == 0.0)
    A=randn(size(A));
  end
  ret = A;
  [M, N]=GramSchmidt(ret, A);
end


