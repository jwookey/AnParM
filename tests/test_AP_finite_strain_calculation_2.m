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
         [FST,c,R,phV(istep)] = update_FST(FST,vgrad,dt,R) ;   
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


function [FST_new,c,R,ph] = update_FST(FST,vgrad,t,PrevR) ;
% update the finite strain tensor based on an applied velocity gradient
% tensor, following McKenzie (1979).
%
% Inputs: 
%    FST    : Previous finite strain tensor (3x3)
%    vgrad  : velocity gradient tensor (3x3)
%    t      : deformation time
%    PrevR  : the previous rotation matrix (3x3)*
%
% Outputs:
%    FST_new : updated finite-strain tensor
%    c       : axis lengths of the (new) finite strain ellipse
%    R       : rotation matrix required to map a vector/tensor in
%              the original cartesian reference frame onto that
%              of the deformed FSE.
%    ph      : minimum angle of this rotation
%
%  *The code calculates the principle axes of the finite strain ellipse, and 
%  calculates the arrangement of these axes which imply the smallest rotation
%  from the previous rotation. This is potentially rather a rate-limiting step.

   
   % form A and B matrix
   A = eye(3,3)-(0.5*t)*vgrad ;
   B = eye(3,3)+(0.5*t)*vgrad ;
   
   % calculate the resulting Finite Strain Tensor.
   FST_new = inv(A)*B*FST ;
   
   % calculate the Cauchy deformation tensor.
   CDT = FST_new'*FST_new ;
   [EIVEC, EIVAL] = eig(CDT) ;
   
   % find the orientation of principle axes which most closely matches
   % the previous axes, i.e., the smallest possible rotation.
   [index]=SortPrincipleAxes(EIVEC,PrevR) 
   
   % index contains minimum distance axes.
   for ii=1:3
      R(:,ii) = sign(index(ii)).*EIVEC(:,abs(index(ii))) ;
      c(ii)=sqrt(EIVAL(abs(index(ii)),abs(index(ii)))) ;
   end
   
   % invert to get the rotation matrix.
   R=inv(R) ;

   ph=-acosd(dot([1 0 0]',R(:,1))) ;
    
end

function [index]=SortPrincipleAxes(E,EPrev) ;
% find the arrangement of the columns of 3x3 matrix E which represents
% the smallest overall rotation from the previous
% negative indices indicate the column values should be reversed.

   index = [0 0 0] ;
   found = [0 0 0] ;

   for ii=1:3
      angle = acosd(dot(E(:,ii),EPrev(:,1))) ;

      if angle<45
         index(1)=ii; break ;
      elseif angle>135
         index(1)=-ii; break ;
      end 
   end

   found(abs(index(1)))=1 ;
   for ii=find(found==0)
      angle = acosd(dot(E(:,ii),EPrev(:,2))) ;
      if angle<45
         index(2)=ii; break ;
      elseif angle>135
         index(2)=-ii; break ;
      end 
   end

   found(abs(index(2)))=1 ;
   for ii=find(found==0)
      angle = acosd(dot(E(:,ii),EPrev(:,3))) ;

      if angle<45
         index(3)=ii; break ;
      elseif angle>135
         index(3)=-ii; break ;
      end 
   end
end  




function [vgrad,ph,r12,r13,r23] = analytical_simple_shear(dt,nstep)
   % from Ribe and Yu, 1991
   
   vgrad = zeros(nstep,3,3) ;
   
   ph = -0.5*atand(([1:nstep]-1)*dt) 
   
   for istep=1:nstep
      E = [ cosd(2*ph(istep))  sind(2*ph(istep)) 0 ; ...
            sind(2*ph(istep)) -cosd(2*ph(istep)) 0 ; ...
            0                 0                0 ] ;
      
      W = [ 0 -1 0 ; 1 0 0 ; 0 0 0 ] ;
      vgrad(istep,:,:) = E + W ;
   end
   
   r12=2*asinh(dt*[0:nstep-1]) ;
   r23=-asinh(dt*[0:nstep-1]) ;
   r13=asinh(dt*[0:nstep-1]) ;
   
end
   