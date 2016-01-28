% Test AnPar method for analytical simple shear example.
function test_AP_finite_strain_calculation_2(nstep)
      
   % power law exponent
   rn = 3.5 ;
   
   % vgrad
   vgrad = [1.0  -1 0 ; 1 -1 0 ; 0 0 0] ;
   %vgrad = zeros(3,3) ;
   
   % CRSS
   tau = [0.3333 0.6667 1.0] ;
   
      
   % time step
   dt = 0.05 ;
   %dt = 0.01 ;
   
   FST = eye(3,3) ; % initial finite strain tensor
   M = eye(3,3) ;
      
   
   for istep = 1:nstep


      
      % rotate the vgrad into the *current* FSE frame: this is breaking it(?)
      vgradR = M*vgrad*M' ;
      
%     ** calc point
      vgradV(istep,:,:) = vgradR ;
      

      
      % calculate the rotation and axes change in FSE
      [FST_new,R,chi] = finitestrain(FST,vgrad,dt) ;
                              
      % plot the finite strain ellipse
      %ivec = R * ivec ;
            
      % accrue the rotation
      M = R*M ;
      
   end
   

   
end


function [FST_new,R,chi] = finitestrain(FST,vgrad,t) ;
   
   % form A and B matrix
   A = eye(3,3)-t*vgrad ;
   B = eye(3,3)+t*vgrad ;
   
   % calculate the resulting Finite Strain Tensor.
   FST_new = inv(A)*B*FST ;
   
   % calculate the Cauchy deformation tensor.
   CDT = inv(FST_new)'*inv(FST_new) 
   
   [EIVEC, EIVAL] = eig(CDT)
   
   % find the arrangement of axes which is closest to the previous axes
   
   c_RAW = [EIVAL(1,1) EIVAL(2,2) EIVAL(3,3)] ;
   [c, IND] = sort(c_RAW,2, 'descend') ;
   R = EIVEC ; % for dimensioning
   for i=1:3
       R(:,i) = EIVEC(:,IND(i)) ;
   end
   
    
end

%function [r12,r23,r13] = finitestrain(vgrad,t) ;
%   
%   % Strain Rate Tensor
%   bige = zeros(3,3) ;
%   for i = 1:3
%      for j = 1:3
%         bige(i,j) = (vgrad(i,j) + vgrad(j,i))/2.d0 ;
%      end
%   end
%   
%   % finite strain ellipse
%   c1 = exp(bige(1,1)*t) ;
%   c2 = exp(bige(2,2)*t) ;
%   c3 = exp(bige(3,3)*t) ;
%      
%   % calculate log ratios
%   r12 = log(c1/c2) ;
%   r23 = log(c2/c3) ;
%   r13 = log(c1/c3) ;   
%   
%   % equivalent strain
%   eijsq = sum(sum(bige.^2)) ;
%   epseq = sqrt(2*eijsq/3)*t ;
%   disp(epseq)
%   
%end






   
   