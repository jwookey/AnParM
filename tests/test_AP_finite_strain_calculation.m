% Test AnPar method for analytical simple shear example.
function test_AP_finite_strain_calculation(nstep)
      
   % power law exponent
   rn = 3.5 ;
   
   % vgrad
   vgrad = [1.0  -1 0 ; 1 -1 0 ; 0 0 0] ;
   % vgrad = [1.0  0 -1 ; 0 0 0 ; 1 0 -1 ] ;

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
   
   [vgradA,phA,r12A,r13A,r23A] = analytical_simple_shear(dt,nstep);
   
   ph = 0 ;   
   
   vgradAR = vgradA ;
   
   FSE = [1.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 1.5];
   FSE = eye(3);
   %FES = [sqrt(2) 0.0 sqrt(2) ;0 1 0; -sqrt(2) 0.0 sqrt(2)]
   det(FSE)
   
   for istep = 1:nstep

      % Update the texture at the current point 
      % for this we need the lengths of the FSE
      % and the rotation onto the reference frame
      % - i.e. the eigenvalues and vectors respctivly
      
      % [c, R] = principAxes(FSE);
      if istep == 1
          R = eye(3);
          c = [1 1 1];
          ph = 0;
      else
          [c, R, ph] = principAxes_XYplaneStrain(FSE);
      end


      % Put the Velocity gradient tensor onto the right 
      % frame of reference.

      vgradR = R'*vgrad*R;
      
      % calculate log ratios
      r12 = log(c(1)/c(2)) ;
      r23 = log(c(2)/c(3)) ;
      r13 = log(c(1)/c(3)) ; 
      
      % Update the texture here...
      
      % [cnew,R,ph] = finitestrain(c,vgradR,dt);
      
      % Update the FSE for this strain incrememnt
      [FSE] = update_finitestrain_McK(FSE,vgrad,dt);
      
      
      % Store some data 
      
      vgradV(istep,:,:) = vgradR ;
      
      r12V(istep) = r12 ;
      r23V(istep) = r23 ;
      r13V(istep) = r13 ;
      
      phV(istep) = ph;
      %if istep==1
      %   phV(istep) = 0 ;
      %else 
      %   phV(istep) = phV(istep-1)+ph ;
      %end
      
      
                            
   end
   

   % now compare
   
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
   
   vgradA(1,:,:)
   vgradV(1,:,:)
   
   vgradA(2,:,:)
   vgradV(2,:,:)
   
   vgradA(3,:,:)
   vgradV(3,:,:)
   
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

function [cnew,R,ph] = finitestrain(c,vgrad,t)
   
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


function [FSE] = update_finitestrain_McK(FSE,vgrad,t) 
   
   % McKenzie '79 eq 8-10
   A = eye(3) - (0.5*t)*vgrad;
   B = eye(3) + (0.5*t)*vgrad;
   FSE = inv(A)*B*FSE;
   
end   

function [c, R] = principAxes(F)
    % Return the lenghts of the principal 
    % axes, c, and the orentations, R, of the 
    % finite strain ellipsoid assoceated 
    % with the deformation gradient tensor, F.
    
    % We need the Eigenvectors and Eigenvalues
    % of the left stretch tensor, V, which is 
    % derived from the Left Cauchy-Green deformation
    % tensor. And we need these to be sorted.

    % Left Cauchy-Green deformation tensor:
    B = F*F';
    
    % Do we need this rotation?
    % Aparantly not...
%     RR = B\sqrtm(F)
%     iRR = inv(RR)
%     tRR = transpose(RR)
    
    % Now, B = V^2, so we can could take the
    % (matrix) square root sqrtm before finding the
    % eigenvalues and vectors, but it is probably 
    % quicker to find the eigenvectors and eigenvalues
    % of B. These are also the eigenvectors of V and the
    % eignevalues of V are just the (element wise) square 
    % root of the eigenvalues of B.
    
    [EIVEC, EIVAL] = eig(B);
    c_RAW = [sqrt(EIVAL(1,1)) sqrt(EIVAL(2,2)) sqrt(EIVAL(3,3))] ;
    [c, IND] = sort(c_RAW,2, 'descend') ;
    R = EIVEC ; % for dimensioning
    for i=1:3
        R(:,i) = EIVEC(:,IND(i)) ;
    end

end

function [c, R, ph] = principAxes_XYplaneStrain(F)

    [c, R] = principAxes(F);
    % For plane strain (e.g. analytical simple shear) we want the 
    % rotation around the X3 axis (i.e. c(3) always = 1)
    % hence swap c(2) and c(3) and the rotation matrix
    c_old3 = c(3);
    c(3) = c(2);
    c(2) = c_old3;
    R_old3 = R(:,3);
    R(:,3) = R(:,2);
    R(:,2) = R_old3;
    
    % Rotation between two frames - plane strain only.
    ph=-acos(-R(1,1));
    %if ph >= pi
    %    ph = ph - pi
    %end
    

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






   
   