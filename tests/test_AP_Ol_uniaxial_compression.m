% Uniaxial compression test of AnParM.
function test_AP_Ol_uniaxial_compression()
   
   % set number of steps
   istepmax = 30 ;
         
   % velocity gradient tensor
   vgrad = [0.5  0.0  0.0 ; ...
            0.0  0.5  0.0 ; ...
            0.0  0.0 -1.0 ] ; % uniaxial compression (x3-axis)
    
   % power law exponent
   rn = 3.5 ;
   
   % CRSS
   tau = [0.3333 0.6667 1.0] ;
      
   % timestep
   dt = 0.02165 ;
      
   % build random texture
   [ texture ] = MVT_make_random_texture( 2000 ) ;
      
   for istep = 1:istepmax
      
      % get FSE axes ratios
      [r12,r23,r13] = finitestrain(vgrad,(istep-1)*dt) ;
      
      % run a texture calculation step
      [new_texture] = ...
            AP_Ol_texture_update(texture,rn,tau,vgrad,r12,r23,r13,dt) ;

      % feed the updated texure back in.                   
      texture = new_texture ;
   end
   
   % output the final texture    
   MVT_write_VPSC_file('uniaxial_compression.out', ...
      texture, 'Uniaxial compression output') ;
      
   % plot with MTEX
   MVT_olivine_pole_from_vpsc('uniaxial_compression.out','scale',[0 17], ...
      'writefile','unaxial_compression','png') ;   

end


function [r12,r23,r13] = finitestrain(vgrad,t) ;
   
   % Strain Rate Tensor
   bige = zeros(3,3) ;
   for i = 1:3
      for j = 1:3
         bige(i,j) = (vgrad(i,j) + vgrad(j,i))/2.d0 ;
      end
   end
   
   % finite strain ellipse
   c1 = exp(bige(1,1)*t) ;
   c2 = exp(bige(2,2)*t) ;
   c3 = exp(bige(3,3)*t) ;
   
   % calculate log ratios
   r12 = log(c1/c2) ;
   r23 = log(c2/c3) ;
   r13 = log(c1/c3) ;   
   
   % equivalent strain
   eijsq = sum(sum(bige.^2)) ;
   epseq = sqrt(2*eijsq/3)*t ;
   
   %disp(epseq)
   
end

   
   