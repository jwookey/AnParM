% AP_velocity_simple_shear - A example function for an arbitary flow field.
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage:
%
%     [V] = AP_velocity_simple_shear(X, t)
%        Input parameters:
%           X : position. 3D position in space.
%           t : time. This example is for time invariant flow, so this 
%               argument has no effect.
%
%        Output parameters:
%           V : 3D velocity at point X and time t.
%
%  A simple analytically tractable example for testing.
%
% Reference.
% ~~~~~~~~~~
%
%  Goulding, NG, Ribe, NM, Castelnau, O., Walker A. and Wookey, J. Analytical
%     parameterization of self-consistent polycrystal mechanics: Fast calculation 
%     of upper mantle anisotropy. Geophys. J. Int., in press.
%
% See also: 

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

function [V] = AP_velocity_simple_shear(X, t)

    assert(all(size(X) == [3 1]), 'AP:argChk', ... 
        'Location must be a 3 vector')
    assert(all(size(t) == [1 1]), 'AP:argChk', ... 
        'Time must be a 3 scalar')
    
    dvy_dx = 1.0;
    V = zeros(3,1);
    V(1) = X(2)*dvy_dx;

end