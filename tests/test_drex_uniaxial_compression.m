% Uniaxial compression test of AnParM.
function test_drex_uniaxial_compression()
   
   % set number of steps
   istepmax = 30 ;
         
   % velocity gradient tensor
   vgrad = [0.5  0.0  0.0 ; ...
            0.0  0.5  0.0 ; ...
            0.0  0.0 -1.0 ] ; % uniaxial compression (x3-axis)
    
   % power law exponent
   rn = 3.5 ;
   
   % CRSS
   tau = [1 2 3 1000] ;
      
   % timestep
   dt = 0.02165 ;
   
   % DRex params
   Mob = 125.0; % Intrinsic grain boundary mobility
   Xol = 1.0; % Volume fraction of olivine
   chi = 0.3; % Threshold for grain boundary sliding
   lambda = 5.0; % Nucleation rate 

   % build random texture
   [ texture ] = MVT_make_random_texture( 2000 ) ;
   % volume fraction for grains
   odf = ones(2000, 1) .* 1/2000;
   rt = zeros(2000, 1);
   
   % FIXME FIXME FIXME: need to calculate a reference strain rate...
   % probably belongsd inside the drex driver... 
   % but I need to think... but for now pick a number so this
   % does something sensible
   epsnot = 0.05;
   
   for istep = 1:istepmax
      
      % run a texture calculation step with DRex
      [new_texture, odf, rt] = ...
            AP_drex_texture_update(texture,rn,tau,vgrad,epsnot, dt, ...
                           Mob, Xol, chi, lambda, odf, rt) ;

      % feed the updated texure back in.                   
      texture = new_texture ;
   end

   
   % FIXME We need a way to store and record the volume fraction
   % in the VPSC file (it's column 4) and account for this in the
   % pole figure (not sure how easy this is)
   
   % rotate the texture FoR to match Goulding et al Figure 7 (vertical compression axis)
   [texture]=AP_rotate_texture_Euler(texture,0,90,90) ;
   
   % output the final texture    
   MVT_write_VPSC_file('uniaxial_compression.out', ...
      texture, 'Uniaxial compression output') ;
      
   % plot with MTEX
   MVT_olivine_pole_from_vpsc('uniaxial_compression.out','scale',[0 17], ...
      'writefile','uniaxial_compression','png') ;   

end
