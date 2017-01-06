% AP_strain_path - Convert a set of velocity gradient tensors and times to 
%                  strain path data compatible with AnPar texture modelling
%
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage:
%
%     [vgradR,R,r12,r23,r13]=AP_strain_path(vgrad,time)
%
%        Input parameters:
%             vgrad : List of velocity gradient tensors (N x 3 x 3). 
%              time : Time associated with each VGT, this should be
%                     the same length as vgrad
%
%        Output parameters:
%            vgradR : (N x 3 x 3) VGT list rotated into the reference frame of
%                     the evolving finite-strain ellipse. 
%                 R : (N x 3 x 3) Rotation matrices which perform these rotations.
%       r12,r23,r13 : (N) Axis ratios of the evolving finite strain ellipse.
%
%
%	The sample is assumed to have no initial deformation. The AnPar method requires
% that the velocity gradient tensor supplied is in the same frame of reference
% as the finite strain ellipse.

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

function [vgradR,R,r12,r23,r13]=AP_strain_path(vgrad,time)
	% do some input checking
	s_vgrad = size(vgrad) ;
	nstep = s_vgrad(1) ;

	if nstep~=length(time)
		error('VGT vector input not same length as time.')
	end

	vgradR = zeros(nstep,3,3) ;
	R = zeros(nstep,3,3) ;
	r12 = zeros(nstep,1) ;
	r23 = zeros(nstep,1) ;
	r13 = zeros(nstep,1) ;

   FST = eye(3,3) ; % initial finite strain tensor (undeformed)

	for istep = 1:nstep
		if istep == 1
          R(istep,:,:) = eye(3);
          c=[1 1 1] ;
    else
      % calculate the rotation and axes change in FST
      [FST,c,R(istep,:,:),~] = AP_update_finite_strain(FST,...
                        squeeze(vgrad(istep,:,:)),...
         								time(istep)-time(istep-1),... % dt
         								squeeze(R(istep-1,:,:)) ) ;   
    end

      % rotate the vgrad into the *current* FSE frame
      vgradR(istep,:,:) = squeeze(    R(istep,:,:))' * ...
                          squeeze(vgrad(istep,:,:))  * ...
                          squeeze(    R(istep,:,:))  ;
      
      % calculate FSE axis length log ratios
      r12(istep) = log(c(1)/c(2)) ;
      r23(istep) = log(c(2)/c(3)) ;
      r13(istep) = log(c(1)/c(3)) ; 

	end

end