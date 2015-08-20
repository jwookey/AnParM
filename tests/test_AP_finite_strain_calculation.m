% Test AnPar method for analytical simple shear example.
function test_AP_Ol_simple_shear_2(nstep)
      
   % power law exponent
   rn = 3.5 ;
   
   % vgrad
   vgrad = [1.0  -1 0 ; 1 -1 0 ; 0 0 0] ;

   % CRSS
   tau = [0.3333 0.6667 1.0] ;
   
      
   % time step
   dt = 0.02165 ;
   
   c = [1 1 1] ; % initial FSE
   M = eye(3,3) ;
      
   ivec = [0 ; 1 ; 0] ;
   
   vgradV = zeros(nstep,3,3) ;
   r12V = zeros(1,nstep) ;
   r13V = zeros(1,nstep) ;
   r23V = zeros(1,nstep) ;
   
   phV = zeros(1,nstep) ;
   
   [vgradA,phA,r12A,r13A,r23A] = analytical_simple_shear(dt,nstep)
   
   ph = 0 ;   
   
   vgradAR = vgradA ;
   
   for istep = 1:nstep

      % RMA = MS_rotM(0,0,phA(istep)*(180/pi)) ; % rotate back to reference frame
      %vgradAR(istep,:,:) = RMA*squeeze(vgradA(istep,:,:))*RMA' ;

      % calculate log ratios
      r12 = log(c(1)/c(2)) ;
      r23 = log(c(2)/c(3)) ;
      r13 = log(c(1)/c(3)) ;      
      
      %AP_plot_FSE(c,M)
      
      % rotate the vgrad into the *current* FSE frame: this is breaking it(?)
      vgradR = M*vgrad*M' ;
      
%     ** calc point
      vgradV(istep,:,:) = vgradR ;
      
      r12V(istep) = r12 ;
      r23V(istep) = r23 ;
      r13V(istep) = r13 ;
      
      if istep==1
         phV(istep) = 0 ;
      else 
         phV(istep) = phV(istep-1)+ph ;
      end
      
      % calculate the rotation and axes change in FSE
      [cnew,R,ph] = finitestrain(c,vgrad,dt) ;
                        
      c=cnew ;
      
      % plot the finite strain ellipse
      %ivec = R * ivec ;
            
      % accrue the rotation
      M = R*M ;
      
   end
   
   % now compare
   
   %vgrad
   figure('Position',[1 1 800 1400]) ;
   subplot(3,1,1)
   title('phi')
   hold on
   plot([1:nstep],phV,'k-') ;
   plot([1:nstep],phA,'k--') ;

   subplot(3,1,2)
   title('vgrad') ;
   hold on
   plot([1:nstep],vgradV(:,1,1),'k-') ;
   plot([1:nstep],vgradA(:,1,1),'k--') ;
   %plot([1:nstep],vgradAR(:,1,1),'k:') ;

   plot([1:nstep],vgradV(:,1,2),'r-') ;
   plot([1:nstep],vgradA(:,1,2),'r--') ;
   %plot([1:nstep],vgradAR(:,1,2),'r:') ;

   plot([1:nstep],vgradV(:,2,1),'g-') ;
   plot([1:nstep],vgradA(:,2,1),'g--') ;
   %plot([1:nstep],vgradAR(:,2,1),'g:') ;

   plot([1:nstep],vgradV(:,2,2),'b-') ;
   plot([1:nstep],vgradA(:,2,2),'b--') ;
   %plot([1:nstep],vgradAR(:,2,2),'b:') ;

   subplot(3,1,3)
   title('Rxx')
   hold on
   plot([1:nstep],r12V,'k-') ;
   plot([1:nstep],r12A,'k--') ;
   plot([1:nstep],r23V,'b-') ;
   plot([1:nstep],r23A,'b--') ;
   plot([1:nstep],r13V,'r-') ;
   plot([1:nstep],r13A,'r--') ;
   
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

function [cnew,R,ph] = finitestrain(c,vgrad,t) ;
   
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
            (vgrad(2,1) - vgrad(1,2))/2] ;
   
   
   
   % form rotation matrix
   [ R ] = MS_rotM( bigom(1)*t*(180/pi), ...
                    bigom(2)*t*(180/pi), ...
                    bigom(3)*t*(180/pi) ) ;
   
   ph=-bigom(3)*t ;
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






   
   