% AP_rotate_texture_Euler - Change the reference frame for a texture represented
%                           by a list of Euler angles into a new reference frame
%                           (itself defined by 3 Euler angles.) 
%
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage:
%
%     [new_texture] = [rotated_texture]=AP_rotate_texture(texture,phi1,Phi,phi2)
%
%        Input parameters:
%           texture : Initial texture to be updated. This is an 3 x M list of
%                     Euler angles in degrees, representing the initial 
%                     orientation of the M crystals.
%     phi1,Phi,phi2 : Euler angles representing the new frame of reference (in
%                     degrees, Bunge notation).
%
%        Output parameters:
%   rotated_texture : the updated texture. An M x 3 list of Euler angles in
%                     degrees, representing the orientation of the 
%                     crystals in the new reference frame.
%
% This uses routines from microtexm (MVT_rot_from_Euler and MVT_rot_to_Euler). 
%
% Notes:
%
%    In the context of this function the term 'Euler angle' is used in 
%    the sense common in texture analysis not computer graphics etc. where
%    the three angles are rotations around a fixes axis system (see 
%    MS_rot3 for this case). In the (Bunge) convention used here the
%    relationship between the orientation of a set of "crystal" axes is 
%    considered with respect to an external "sample" axis system. This 
%    is described by considering a moving reference frame which starts out 
%    parallel to the sample axes. This moving reference frame is first 
%    rotated anticlockwise about its z axis by phi1 degrees. A second 
%    anticlockwise rotation of Phi degrees is made about the frame's x axis 
%    (which, at this point, is no longer parallel to the sample x axis).
%    Finally, a second anticlockwise rotation of phi2 degrees about the 
%    frame's z axis is made to bring the moving reference into line with
%    the crystal axes. See Figure 5 of Bunge 1985 for a graphical example.
%

% Copyright (c) 2015, Neil Goulding, Neil Ribe, Olivier Castelnau, 
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

function [rotated_texture]=AP_rotate_texture(texture,phi1,Phi,phi2)
   
   % Form a rotation matrix describing the rotation of
   % the texture
   g_2 = MVT_rot_from_Euler(phi1,Phi,phi2);
   
   % apply it to the crystals
   rotated_texture = zeros(size(texture));
   
   for i=1:length(texture) ;
      g_1 = MVT_rot_from_Euler(texture(1,i), texture(2,i), texture(3,i)) ; 
      g_1 = g_1*g_2' ;
      [phi1, Phi, phi2] = MVT_rot_to_Euler(g_1) ;
      rotated_texture(1,i) = phi1;
      rotated_texture(2,i) = Phi;
      rotated_texture(3,i) = phi2;         
   end
   
end