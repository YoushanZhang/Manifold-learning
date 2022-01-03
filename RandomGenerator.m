function data=  RandomGenerator(mean, sigma, method) 

if nargin<3
    method = 'Sphere';
end

randomList=zeros(size(mean));
 for i=1:1:3
         r2 = 0;
      while(r2 > 1.0 | r2 == 0)
     
        x = randn;
        y = randn;
        r2 = x * x + y * y;
      end
     randomNumber = sigma * y * sqrt(-2.0 * log(r2) / r2);
     randomList(i) = randomNumber;
 end
 temp=listToTangent(randomList,method);
 randomTangent=projectTangent( mean, temp,method);
 data = exp_Map(mean, randomTangent,method);
    
      
