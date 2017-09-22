% AP_drex_texture_update - Update an input  texture (as a set of Euler 
%                          angles) based on a velocity gradient tensor, 
%                          and a set of slip system parameters using the
%                          'DRex' approach.
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //     
%
% FIXME: Usage and documentation here!

% Copyright (c) 2016, Neil Goulding, Neil Ribe, Olivier Castelnau, 
%                     Andrew Walker and James Wookey
%
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in the 
%      documentation and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names of its 
%      contributors may be used to endorse or promote products derived from 
%      this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% FIXME: this assumes we have four slip systems. Somewhere I have a 
% version where this is a loop ovrt the size (it's an OpenMP version)
% need to find this!

function [ new_texture, odf, rt ] = ...
   AP_drex_texture_update( texture, rn, tau, vgrad, epsnot, dt, ...
                           Mob, Xol, chi, lambda, odfi, rt)

   % NB - different frame of reference to ANPAR - external frame not FSE
   % frame!
   
   % FIXME: what about macroscopic rotation??

   ngrains = length(texture) ;
   
   % Should we pass odf and update at the end for each step??
   assert(length(odfi) == ngrains);
   assert(length(rt) == ngrains);

   % Strain Rate Tensor
   bige = zeros(3,3) ;
   for i = 1:3
      for j = 1:3
         bige(i,j) = (vgrad(i,j) + vgrad(j,i))/2.d0 ;
      end
   end 

   % Setup array of direction cosines
   acsi = zeros(ngrains,3,3);
   for i = 1:ngrains
       acsi(i,:,:) = MVT_rot_from_Euler(...
           texture(1,i), texture(2,i), texture(3,i));
   end
   
   % Other initalisation
   dotacs = zeros(ngrains,3,3);
   dotodf = zeros(ngrains);
   new_texture = zeros(3, ngrains);
   odf = zeros(ngrains,1);
   
   % Called 'alt' in DRex - and shared with AnPar code too
   % FIXME: cleanup needed.
   eps = zeros([3 3 3]) ;
   eps(1,2,3) = 1 ;
   eps(2,3,1) = 1 ;
   eps(3,1,2) = 1 ;
   eps(2,1,3) = -1 ;
   eps(3,2,1) = -1 ;
   eps(1,3,2) = -1 ;
   
   
   
   % Stuff below this is from 'deriv.f'
   
   
   % dimensionless strainrate and velgrad tensors
   lx = vgrad / epsnot;
   ex = bige / epsnot;
   
   ngrains = length(texture) ;
   
   % Plastic deformation and dynamic recrystallization
   for i = 1:ngrains
   
       bigi = zeros(4,1);
       gam = zeros(4,1);
       g = zeros(3,3);
       
       % Big I in eq 5 of Kaminski and Ribe
       % note implicit sum on i and j
       % ! NB Relationship between slip systems and crystal axis is
       % ! implicit here. e.g. for system #1 assumes that the slip
       % ! plane normal is parallel to the b axis and the slip direction
       % ! is parallel to the a axis. This is right for (010)[100] in a
       % ! orthorhombic crystal, but will need generalizing for the general
       % ! case. See DRex Note A in my notebook.
       for i1 = 1:3
           for i2 = 1:3
               bigi(1) = bigi(1)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,2,i2);
               bigi(2) = bigi(2)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,3,i2);
               bigi(3) = bigi(3)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,2,i2);
               bigi(4) = bigi(4)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,1,i2);
           end
       end
       
       % Quotients
       q = bigi./tau';
   
       % Sort: qab(imax) > qab(iint) > qab(imin) > qab(iinac)
       qab = abs(q);
       [~, ti] = sort(qab);
       imax = ti(4);
       iint = ti(3);
       imin = ti(2);
       iinac = ti(1);
       
       % Calculate weighting factors
       % This is eq. 5 of Kaminski and Ribe 2001
       gam(imax) = 1.0;
       
       rat = tau(imax)/bigi(imax);
       qint = rat*bigi(iint)/tau(iint);
       qmin = rat*bigi(imin)/tau(imin);
       sn1 = rn - 1.0; % we call stressexp rn!
       
       gam(iint) = qint*(abs(qint))^sn1;
       gam(imin) = qmin*(abs(qmin))^sn1;
       gam(iinac) = 0.0;
       
       % Calculate G tensor
       % This is e.q. 4 of Kaminski and Ribe 2001
       % See note avove re relationship between slip system. gam*s( is beta
       % ub tge equation
       for i1 = 1:3
           for i2 = 1:3
               g(i1, i2) = 2.0*(gam(1)*acsi(i,1,i1)*acsi(i,2,i2) + ...
                                gam(2)*acsi(i,1,i1)*acsi(i,3,i2) + ...
                                gam(3)*acsi(i,3,i1)*acsi(i,2,i2) + ...
                                gam(4)*acsi(i,3,i1)*acsi(i,1,i2));
           end
       end
       
       % Calculate strain rate on the softest slip system
       % This is e.q. 7 of Kaminski and Ribe 2001 but j is 'i', k is 'j'
       % and i2 is 'l' (note the cyclic thing). Summing repeated terms. 
       % gam0 is a scalar for each crystal...
       R1 = 0.0;
       R2 = 0.0;
       for j = 1:3
           i2 = j + 2;
           if (i2 > 3) 
               i2 = i2 - 3;
           end
           R1 = R1 + (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j));
           R2 = R2 + (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j));
           for k = 1:3
               R1 = R1 + 2.0*g(j,k)*g(j,k);
               R2 = R2 + 2.0*lx(j,k)*g(k,k);
           end
       end
       gam0 = R2/R1;
       
       % Calculate the rotation rate (vector)
       % This is e.q. 8 of Kaminski and Ribe 2001
       rot(3) = (lx(2,1)-lx(1,2))/2.0 - (g(2,1)-g(1,2))/2.0 * gam0;
       rot(2) = (lx(1,3)-lx(3,1))/2.0 - (g(1,3)-g(3,1))/2.0 * gam0;
       rot(1) = (lx(3,2)-lx(2,3))/2.0 - (g(3,2)-g(2,3))/2.0 * gam0;

       % Derivative of the matrix of direction cosine (e.q. 9 of Kamiski 
       % and Ribe 2001
       
       % Fixme - need dotacs to be zeros like acs...
       for i1 = 1:3
           for i2 = 1:3
               for i3 = 1:3
                   for i4 = 1:3
                       dotacs(i,i1,i2) = dotacs(i,i1,i2) + ...
                           eps(i2,i3,i4)*acsi(i,i1,i4)*rot(i3);
                   end
               end
           end
       end
       
       % Dislocation density calculation
       % e.q. 10 and 11 of Kaminski, Ribe and Browaeys
       rt1 = tau(1)^(1.5-rn) * abs(gam(1)*gam0)^(1.5/rn);
       rt2 = tau(2)^(1.5-rn) * abs(gam(2)*gam0)^(1.5/rn);
       rt3 = tau(3)^(1.5-rn) * abs(gam(3)*gam0)^(1.5/rn);
       
       rt(i) = rt1*exp(-lambda*rt1^2) + ...
               rt2*exp(-lambda*rt2^2) + ...
               rt3*exp(-lambda*rt3^2);
           
       % Grain boundary sliding for small grains
       % e.q. 15 and 16 of Kaminski, Ribe and Browaeys
       if (odfi(i) < chi/ngrains)
           dotacs(i, :, :) = 0;
           rt(i) = 0;
       end
       
   end % of loop over all grains
   
   % Volume av energy (e.q. 15 and 16 of Kaminski, Ribe and Brownaeys
   emean = sum(odfi.*rt);
   
   % Update volume fraction by grain boundary migtation
   % e.g. 12 of Kaminski, Ribe and Brownaeys
   for i = 1:ngrains
       dotodf(i) = Xol * Mob * odfi(i) * (emean-rt(i));
   end
   
   % Don't do enstatite yet (should be able to turn the above into 
   % a function that works for both phases
   
   
   % Now update the texture. NB: this is done by 4th order RK int
   % in 'strain.f' in DRex, here we just take a step...
   % FIXME: is this OK? What about the strain rate normalisation
   % above. Do we need the *dt?? 
   % NB: rt is a value and not a derivative
   g = zeros(3,3);
   for i = 1:ngrains
       acsi(i,:,:) = acsi(i,:,:) + dt.*dotacs(i,:,:);
       g(:,:) = acsi(i,:,:);
       [phi1, Phi, phi2] = MVT_rot_to_Euler(g);
       new_texture(1,i) = phi1;
       new_texture(2,i) = Phi;
       new_texture(3,i) = phi2;  
       odf(i) = odfi(i) + dt.*dotodf(i);
   end
   

end


