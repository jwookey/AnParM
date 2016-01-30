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
   
   % finite strain tensor 
   FST = eye(3,3) ;

   % build random texture
   [ texture ] = MVT_make_random_texture( 2000 ) ;
      
   for istep = 1:istepmax
      
      % get FSE axes ratios
      [FST,r12,r23,r13] = finitestrain(FST,vgrad,dt) ;
      
      % run a texture calculation step
      [new_texture] = ...
            AP_Ol_texture_update(texture,rn,tau,vgrad,r12,r23,r13,dt) ;

      % feed the updated texure back in.                   
      texture = new_texture ;
   end
   
   % rotate the texture FoR to match Goulding et al Figure 7 (vertical compression axis)
   [texture]=AP_rotate_texture_Euler(texture,0,90,90) ;
   
   % output the final texture    
   MVT_write_VPSC_file('uniaxial_compression.out', ...
      texture, 'Uniaxial compression output') ;
      
   % plot with MTEX
   MVT_olivine_pole_from_vpsc('uniaxial_compression.out','scale',[0 17], ...
      'writefile','uniaxial_compression','png') ;   

end


function [FST_new,r12,r23,r13] = finitestrain(FST,vgrad,t) ;
   
   [FST_new,c,~,~] = AP_update_finite_strain(FST,vgrad,t,eye(3,3)) ;

   % calculate log ratios
   r12 = log(c(1)/c(2)) ;
   r23 = log(c(2)/c(3)) ;
   r13 = log(c(1)/c(3)) ; 
   
end

   
   