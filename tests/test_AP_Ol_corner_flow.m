% Try to reproduce figure 10 of Goulding

function test_AP_Ol_corner_flow()

  % Build the strain path using corner flow velocity field
  
  % time vector
  time = [0:0.05:9*0.05 ] ;
  
  % velocity field. Use alpha = 60 degrees (pg 342) and set U_0 to 1.0
  % I suspect U_0 ought to depend on the time step...
  [corner_flow_func] = AP_make_corner_flow(60, 1.0, 'method', 'goulding');

  % FIXME: need to calculate this
  end_location = [0.2 0.2 0.2]';
  
  % Find path through velocity field. We should probably go backwards 
  % using the 'backwards' option here.
  [Xs, Ls] = AP_path_integrate(end_location, 9*0.05, 0.0, -0.05, 0.001,...
      corner_flow_func);
  
  % calculate strain along path as required for AnPar
  [vgradR,R,r12,r23,r13]=AP_strain_path(Ls,time) ;

  vgradR
  
  % power law exponent (but this is 3 for the flow model!)
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
   

  % FIXME: do we need this?
  % rotate the texture FoR to match Goulding et al Figure 8 (horizontal shear plane)
  % This does not seem to be quite right.
  %[texture]=AP_rotate_texture_Euler(texture,0,0,-45) ;
   
  % output the final texture    
  MVT_write_VPSC_file('corner_flow.out', ...
     texture, 'Corner flow output') ;
      
  % plot with MTEX
  MVT_olivine_pole_from_vpsc('corner_flow.out','scale',[0 7], ...
     'writefile','corner_flow','png') 

end