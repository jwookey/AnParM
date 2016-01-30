% Test AnPar method for analytical simple shear example.
function test_AP_Ol_simple_shear_1(nstep)
      
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
   
   [vgradA,phA,r12A,r13A,r23A] = analytical_simple_shear(dt,nstep);

   % build random texture
   [ texture ] = MVT_make_random_texture( 2000 ) ;
   
   ph = 0 ;
   
   for istep = 1:nstep
   
      if istep == 1
          R = eye(3);
          phV(1) = 0;
          c=[1 1 1] ;
      else
         % calculate the rotation and axes change in FST
         [FST,c,R,phV(istep)] = AP_update_finite_strain(FST,vgrad,dt,R) ;   
      end

      % rotate the vgrad into the *current* FSE frame
      vgradR = R'*vgrad*R ;
      
      % calculate log ratios
      r12 = log(c(1)/c(2)) ;
      r23 = log(c(2)/c(3)) ;
      r13 = log(c(1)/c(3)) ; 

      % update the texture
      [new_texture] = ...
            AP_Ol_texture_update(texture,rn,tau,vgradR,r12,r23,r13,dt) ;

      % feed the updated texure back in                   
      texture = new_texture ;

      % save the values at this point
      vgradV(istep,:,:) = vgradR ;
      
      r12V(istep) = r12 ;
      r23V(istep) = r23 ;
      r13V(istep) = r13 ;

      % done
      

   end

   % output the final texture
   % MVT_write_VPSC_file('simple_shear2.out', texture, 'Simple shear output')
   
   % rotate the texture FoR to match Goulding et al Figure 8 (horizontal shear plane)
   % This does not seem to be quite right.
   %[texture]=AP_rotate_texture_Euler(texture,0,0,(ph/pi)*180) ;
   
   % output the final texture    
   %MVT_write_VPSC_file('simple_shear.out', ...
   %   texture, 'Simple shear output') ;
      
   % plot with MTEX
   %MVT_olivine_pole_from_vpsc('simple_shear.out','scale',[0 7], ...
   %   'writefile','simple_shear','png') 

   % compare with the analytical case

   %vgrad
   figure('Position',[1 1 800 1400]) ;
   subplot(3,1,1)
   title('Frame rotation (phi)')
   hold on
   plot([1:nstep],phV,'k-') ;
   plot([1:nstep],phA,'ko') ;
   plot([1:nstep],(-pi/4.0)*ones(nstep, 1), 'k--')
   legend('numerical', 'analytical', 'limit')

   subplot(3,1,2)
   title('vgrad') ;
   hold on
   plot([1:nstep],vgradV(:,1,1),'k-') ;
   plot([1:nstep],vgradA(:,1,1),'ko') ;
   %plot([1:nstep],vgradAR(:,1,1),'k:') ;

   plot([1:nstep],vgradV(:,1,2),'r-') ;
   plot([1:nstep],vgradA(:,1,2),'ro') ;
   %plot([1:nstep],vgradAR(:,1,2),'r:') ;

   plot([1:nstep],vgradV(:,2,1),'g-') ;
   plot([1:nstep],vgradA(:,2,1),'go') ;
   %plot([1:nstep],vgradAR(:,2,1),'g:') ;

   plot([1:nstep],vgradV(:,2,2),'b-') ;
   plot([1:nstep],vgradA(:,2,2),'bo') ;
   %plot([1:nstep],vgradAR(:,2,2),'b:') ;
   legend('D(1,1) - numerical','D(1,1) - analytical', ...
          'D(1,2) - numerical','D(1,2) - analytical', ...
          'D(2,1) - numerical','D(2,1) - analytical', ...
          'D(2,2) - numerical','D(2,2) - analytical')

   subplot(3,1,3)
   title('ln FSE lengths')
   hold on
   plot([1:nstep],r12V,'k-') ;
   plot([1:nstep],r12A,'ko') ;
   plot([1:nstep],r23V,'b-') ;
   plot([1:nstep],r23A,'bo') ;
   plot([1:nstep],r13V,'r-') ;
   plot([1:nstep],r13A,'ro') ;
   legend('r12 - numerical','r12 - analytical', ...
          'r23 - numerical','r23 - analytical', ...
          'r13 - numerical','r13 - analytical')


end         

function [vgrad,ph,r12,r13,r23] = analytical_simple_shear(dt,nstep)
   % from Ribe and Yu, 1991
   
   vgrad = zeros(nstep,3,3) ;
   
   ph = -0.5*atan(([1:nstep]-1)*dt) ;
   
   for istep=1:nstep
      E = [ cos(2*ph(istep))  sin(2*ph(istep)) 0 ; ...
            sin(2*ph(istep)) -cos(2*ph(istep)) 0 ; ...
            0                 0                0 ] ;
      
      W = [ 0 -1 0 ; 1 0 0 ; 0 0 0 ] ;
      vgrad(istep,:,:) = E + W ;
   end
   
   r12=2*asinh(dt*[0:nstep-1]) ;
   r23=-asinh(dt*[0:nstep-1]) ;
   r13=asinh(dt*[0:nstep-1]) ;
   
end

function [vgrad,r12,r23,r13] = setup_VGRAD_FSE(i)
% lookup table for vgrad, FSE from Neil G's Maple code.
% THIS IS NO LONGER USED; see analytical_simple_shear above

%      V11          V12          V21          V22
VXX = [1.0000000000 -1.000000000 1.0000000000 -1.0000000000 ; ...
       0.9997657211 -1.021644928 0.9783550721 -0.9997657211 ; ...
       0.9990638711 -1.043259466 0.9567405344 -0.9990638711 ; ...
       0.9978973988 -1.064813436 0.9351865640 -0.9978973988 ; ...
       0.9962711803 -1.086277084 0.9137229158 -0.9962711803 ; ...
       0.9941919634 -1.107621280 0.8923787200 -0.9941919634 ; ...
       0.9916682900 -1.128817711 0.8711822891 -0.9916682900 ; ...
       0.9887104003 -1.149839061 0.8501609388 -0.9887104003 ; ...
       0.9853301200 -1.170659177 0.8293408232 -0.9853301200 ; ...
       0.9815407322 -1.191253212 0.8087467883 -0.9815407322 ; ...
       0.9773568386 -1.211597756 0.7884022445 -0.9773568386 ; ...
       0.9727942100 -1.231670941 0.7683290588 -0.9727942100 ; ...
       0.9678696322 -1.251452530 0.7485474697 -0.9678696322 ; ...
       0.9626007464 -1.270923980 0.7290760199 -0.9626007464 ; ...
       0.9570058902 -1.290068485 0.7099315147 -0.9570058902 ; ...
       0.9511039390 -1.308871004 0.6911289958 -0.9511039390 ; ...
       0.9449141523 -1.327318262 0.6726817377 -0.9449141523 ; ...
       0.9384560248 -1.345398740 0.6546012601 -0.9384560248 ; ...
       0.9317491461 -1.363102642 0.6368973579 -0.9317491461 ; ...
       0.9248130686 -1.380421856 0.6195781444 -0.9248130686 ; ...
       0.9176671857 -1.397349892 0.6026501085 -0.9176671857 ; ...
       0.9103306220 -1.413881817 0.5861181827 -0.9103306220 ; ...
       0.9028221328 -1.430014182 0.5699858182 -0.9028221328 ; ...
       0.8951600174 -1.445744931 0.5542550693 -0.8951600174 ; ...
       0.8873620432 -1.461073318 0.5389266823 -0.8873620432 ] ;

% construct velocity gradients tensor
vgrad = [VXX(i,1) VXX(i,2) 0.000000 ; ...
         VXX(i,3) VXX(i,4) 0.000000 ;  ...
         0.000000 0.000000 0.000000] ;

%      R12             R23               R13
RXX = [0.00000000000   -0.00000000000    0.00000000000 ; ...
       0.04329661810   -0.02164830905    0.02164830905 ; ...
       0.08657296190   -0.04328648095    0.04328648095 ; ...
       0.12980884240   -0.06490442118    0.06490442118 ; ...
       0.17298424010   -0.08649212003    0.08649212003 ; ...
       0.21607938740   -0.10803969370    0.10803969370 ; ...
       0.25907484820   -0.12953742410    0.12953742410 ; ...
       0.30195159340   -0.15097579670    0.15097579670 ; ...
       0.34469107180   -0.17234553590    0.17234553590 ; ...
       0.38727527640   -0.19363763820    0.19363763820 ; ...
       0.42968680460   -0.21484340230    0.21484340230 ; ...
       0.47190891160   -0.23595445580    0.23595445580 ; ...
       0.51392555860   -0.25696277930    0.25696277930 ; ...
       0.55572145200   -0.27786072600    0.27786072600 ; ...
       0.59728207820   -0.29864103910    0.29864103910 ; ...
       0.63859372940   -0.31929686470    0.31929686470 ; ...
       0.67964352480   -0.33982176240    0.33982176240 ; ...
       0.72041942200   -0.36020971100    0.36020971100 ; ...
       0.76091022620   -0.38045511310    0.38045511310 ; ...
       0.80110559000   -0.40055279500    0.40055279500 ; ...
       0.84099600980   -0.42049800490    0.42049800490 ; ...
       0.88057281640   -0.44028640820    0.44028640820 ; ...
       0.91982816040   -0.45991408020    0.45991408020 ; ...
       0.95875499600   -0.47937749800    0.47937749800 ; ...
       0.99734705740   -0.49867352870    0.49867352870 ] ;

% Finite strain ellipse ratios
r12 = RXX(i,1) ;
r23 = RXX(i,2) ;
r13 = RXX(i,3) ;


end   
   
   
   
   



   
   