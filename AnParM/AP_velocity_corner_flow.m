% AP_velocity_corner_flow - A example function for an arbitary flow field.
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage:
%
%     [V] = AP_velocity_corner_flow(X, t)
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

function [V] = AP_velocity_corner_flow(X, t)

    assert(all(size(X) == [3 1]), 'AP:argChk', ... 
        'Location must be a 3 vector')
    assert(all(size(t) == [1 1]), 'AP:argChk', ... 
        'Time must be a 3 scalar')
    
    theta_0 = 0.5*pi;
    U = 1.0;
    
    % This should be passed in , or it will be slow!
    persistent fdphi_dr;
    persistent fdphi_dtheta;
    if isempty(fdphi_dr)
        [~, fdphi_dr, fdphi_dtheta] = ...
            stream_function_and_derivs_batchelor;
    end
    
    % Convert x, y, z to r, theta, z
    r = sqrt(X(1)^2 + X(2)^2);
    theta = atan2(X(2), X(1));
    
    
    % Calculate velocities from stream function
    % c.f. Batchelor EQ 2.2.10 etc.
    u_r = 1/r * fdphi_dtheta(r, theta, U, theta_0);
    u_theta = -fdphi_dr(r, theta, U, theta_0);
    
    % Convert velocities to cartesian frame
    V = zeros(3,1);
    V(1) = u_r * cos(theta) - u_theta * sin(theta);
    V(2) = u_r * sin(theta) + u_theta * cos(theta);

end

function [fphi, fdphi_dr, fdphi_dtheta] = ...
    stream_function_and_derivs_batchelor
    % Returns Matlab functions for the stream function and the derivatives
    % in phi and r for corner flow following the development of Batchelor
    % (1967) pages 224 - 226.

    theta = sym('theta');
    r = sym('r');
    theta_0 = sym('theta_0');
    U = sym('U');
    
    % EQ 4.8.26
    A = -theta_0^2 * (U/(theta_0^2 - sin(theta_0)^2));
    B = 0; % But we may need this for other boundary conditions
    C = (theta_0 - sin(theta_0) * cos(theta_0)) * ...
        (U/(theta_0^2 - sin(theta_0)^2));
    D = sin(theta_0)^2 * (U/(theta_0^2 - sin(theta_0)^2));
    
    % EQ 4.8.25
    f = A*sin(theta) + B*cos(theta) + C*theta*sin(theta) + ...
        D*theta*cos(theta);
    % EQ 4.8.24
    phi = r * f;
    
    dphi_dr = diff(phi, r);
    dphi_dtheta = diff(phi, theta);
    
    % latex(dphi_dr)
    % latex(dphi_dtheta)
    
    fphi = matlabFunction(phi, 'vars', [r, theta, U, theta_0]);
    fdphi_dr = matlabFunction(dphi_dr, 'vars', [r, theta, U, theta_0]);
    fdphi_dtheta = matlabFunction(dphi_dtheta, ...
        'vars', [r, theta, U, theta_0]);
    
end