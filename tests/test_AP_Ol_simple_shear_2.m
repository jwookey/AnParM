% Test AnPar method for analytical simple shear example.
function AP_odfcalc_Ol_simple_shear_test(istepmax)
      
   % power law exponent
   rn = 3.5 ;
   
   % vgrad
   vgrad = [1 -1 0 ; 1 -1 0 ; 0 0 0] ;
   
   % CRSS
   tau = [0.3333 0.6667 1.0] ;
   
   % number of grains
   ngrains = 2000 ;
      
   % time step
   dt = 0.02165 ;
   
   % number of slip systems
   nslip = 3 ;
   
   % build initial random texture
   [ eulers ] = MVT_make_random_texture( ngrains ) ;
   
   % make it an n x 3 list, and convert to radians
   initial_texture = (pi/180) * eulers' ;
   
   texture = initial_texture ;
   
   c = [1 1 1] ; % initial FSE
   M = eye(3,3) ;
   
   if istepmax>25, error('istepmax must be <=25'), end ;
   
   for istep = 1:istepmax
      % get vgrad and FSE
      %[vgrad,r12,r23,r13] = setup_VGRAD_FSE(istep) ;
      
      % calculate log ratios
      r12 = log(c(1)/c(2)) ;
      r23 = log(c(2)/c(3)) ;
      r13 = log(c(1)/c(3)) ;      
      
      % rotate the vgrad into the *current* FSE frame: this is breaking it.
      vgradR = M*vgrad 
      
      % run a texture calculation step
      [new_texture] = AP_odfcalc_Ol_step(texture,vgrad,...
                      r12,r23,r13,rn,tau,ngrains,dt,nslip) ;

      % feed the updated texure back in                   
      texture = new_texture ;
      
      % calculate the rotation and change in FSE
      [cnew,R] = finitestrain(c,vgrad,dt) ;
            
      c=cnew ;
      
      % accrue the rotation
      M = M*R ;
      
   end
      
   % output the final texture
   texture_out = 180/pi * new_texture' ;
   MVT_write_VPSC_file('simple_shear2.out', texture_out, 'Simple shear output')
   
         
end

function [cnew,R] = finitestrain(c,vgrad,t) ;
   
   % Strain Rate Tensor
   bige = zeros(3,3) ;
   for i = 1:3
      for j = 1:3
         bige(i,j) = (vgrad(i,j) + vgrad(j,i))/2.d0 ;
      end
   end
   
   % finite strain ellipse axes
   cnew = c .* [exp(bige(1,1)*t) exp(bige(2,2)*t) exp(bige(3,3)*t) ] ;
   
   % rotation rates vector
   bigom = [(vgrad(3,2) - vgrad(2,3))/2 ...
            (vgrad(1,3) - vgrad(3,1))/2 ...
            (vgrad(2,1) - vgrad(1,2))/2] ;
   
   % rotation matrix
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






   
   