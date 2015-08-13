% Test AnPar method for analytical simple shear example.
function test_AP_Ol_simple_shear_2(istepmax)
      
   % power law exponent
   rn = 3.5 ;
   
   % vgrad
   vgrad = [1.0  -1 0 ; 1 -1 0 ; 0 0 0] ;
   
   % CRSS
   tau = [0.3333 0.6667 1.0] ;
   
   % number of grains
   ngrains = 2000 ;
      
   % time step
   dt = 0.02165 ;
   
   % build initial random texture
   [ texture ] = MVT_make_random_texture( ngrains ) ;
   
   c = [1 1 1] ; % initial FSE
   M = eye(3,3) ;
   
   if istepmax>25, error('istepmax must be <=25'), end ;
   
   ivec = [0 ; 1 ; 0] ;
   
   for istep = 1:istepmax
      
      % calculate log ratios
      r12 = log(c(1)/c(2)) ;
      r23 = log(c(2)/c(3)) ;
      r13 = log(c(1)/c(3)) ;      
      
      
      % rotate the vgrad into the *current* FSE frame: this is breaking it.
      vgradR = M*vgrad*M' 
      
      fprintf('r12,r23,r13: %8.4f %8.4f %8.4f\n',r12,r23,r13)
      disp('vgradR:')
      disp(vgradR) ;
      
      % run a texture calculation step
      [new_texture] = ...
            AP_Ol_texture_update(texture,rn,tau,vgradR,r12,r23,r13,dt) ;

      % feed the updated texure back in                   
      texture = new_texture ;
      
      % calculate the rotation and change in FSE
      [cnew,R] = finitestrain(c,vgrad,dt) ;
            
      c=cnew ;
      
      % plot the finite strain ellipse
      
      
      ivec = R * ivec ;
            
      % accrue the rotation
      M = R*M ;
      
   end
      
   % output the final texture
   MVT_write_VPSC_file('simple_shear2.out', texture, 'Simple shear output')
   
         
end

function [cnew,R] = finitestrain(c,vgrad,t) ;
   
   % Strain Rate Tensor
   bige = zeros(3,3) ;
   for i = 1:3
      for j = 1:3
         bige(i,j) = (vgrad(i,j) + vgrad(j,i))/2.d0 ;
      end
   end
   
   % finite strain ellipse axes perturbation
   cnew = c .* [exp(bige(1,1)*t) exp(bige(2,2)*t) exp(bige(3,3)*t) ] ;
   

   % rotation rates vector
   bigom = [(vgrad(3,2) - vgrad(2,3))/2 ...
            (vgrad(1,3) - vgrad(3,1))/2 ...
            (vgrad(2,1) - vgrad(1,2))/2] 
   
   % form rotation matrix
   [ R ] = MS_rotM( bigom(1)*t*(180/pi), ...
                    bigom(2)*t*(180/pi), ...
                    bigom(3)*t*(180/pi) ) ;
   
   
   % equivalent strain
   % eijsq = sum(sum(bige.^2)) ;
   % epseq = sqrt(2*eijsq/3)*t ;
   % disp(epseq)
   
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






   
   