% Test AnPar method for analytical simple shear example. 
% This applies the same deformation as test_AP_Ol_simple_shear_1.m, except
% it demonstrates the use of AP_strain_path to create the set of VGTs,
% and demonstrates how to apply AnPar texture updates using them.
%
function test_AP_Ol_simple_shear_2()
   
  % build the strain path, constant (simple shear) velocity gradient tensor
 
  % time vector
  time = [0:0.05:9*0.05 ] ;

  length(time)

  % velocity gradient tensor 
  vgrad = zeros(length(time),3,3) ;
  n=length(time) ;
  for i=1:length(time)
      vgrad(i,:,:) = [1.0  -1 0 ; 1 -1 0 ; 0 0 0] ;
  end

  % calculate strain path details required for AnPar
  [vgradR,R,r12,r23,r13]=AP_strain_path(vgrad,time) ;

  % run the texture update 
  calc_simple_shear(time,vgradR,R,r12,r23,r13)

end

function calc_simple_shear(time,vgradR,R,r12,r23,r13)

   % power law exponent
   rn = 3.5 ;
   
   % CRSS
   tau = [0.3333 0.6667 1.0] ;
   
   % build initial (random) texture
   [ texture ] = MVT_make_random_texture( 2000 ) ;
      
   for i=2:length(time)

      % update the texture
      [new_texture] = ...
            AP_Ol_texture_update(texture,rn,tau,...
                                 squeeze(vgradR(i,:,:)),...
                                 r12(i),r23(i),r13(i),(time(i)-time(i-1))) ;

      % feed the updated texure back in                   
      texture = new_texture ;

   end

   % output the final texture
   MVT_write_VPSC_file('simple_shear2.out', texture, 'Simple shear output')
   
   % rotate the texture FoR to match Goulding et al Figure 8 (horizontal shear plane)
   % This does not seem to be quite right.
   [texture]=AP_rotate_texture_Euler(texture,0,0,-45) ;
   
   % output the final texture    
   MVT_write_VPSC_file('simple_shear.out', ...
      texture, 'Simple shear output') ;
      
   % plot with MTEX
   MVT_olivine_pole_from_vpsc('simple_shear.out','scale',[0 7], ...
      'writefile','simple_shear','png') 


end         




   
   
   
   



   
   