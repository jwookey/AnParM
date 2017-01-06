% AP_path_integrate - Integrate forwards or backwards along a pathline
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage:
%
%     [Xs, Ls] = AP_path_integrate( X0, tstart, tend, dt, dx, @velocity_f)
%        Input parameters:
%                X0 : Initial location.
%            tstart : Starting time.
%              tend : End time
%                dt : Time step for path integration
%                dx : Finite displacment to calculate velocity gradient
%          velocity : function handle allowing evaluation of the velocity
%                     field at a point and time
%
%        Output parameters:
%                Xs : Array of positions of size (3,n).
%                Ls : Array of velocity gradient tensors, size (3,3,n)
%                ts : Array of times, size (1, n)
% 
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

function [Xs, Ls, ts] = AP_path_integrate( X0, tstart, tend, ...
    dt, dx, velocity_f, varargin)

    % Set up defaults
    backwards = 0;
    integrator = @AP_rk4_step;
    
    % Process the optional arguments
    iarg = 1 ;
    while iarg<=(length(varargin))
       switch varargin{iarg}
          case 'backwards'
             backwards = 1 ;
             iarg = iarg + 1 ;
          case 'integrator'
             integrator = varargin{iarg+1} ;
             iarg = iarg + 2 ;
          otherwise 
             error('Unknown flag') ;   
       end   
    end 


    if (backwards)
        
        % Get points going backwards by using a modified 
        % rk4 step, inverting the end time and dt.
        [Xs, Ls_f, ts] = AP_path_integrate( X0, tstart, -tend, ...
            -dt, dx, velocity_f, 'integrator', @AP_rk4_step_backwards);
        
        ts = fliplr(ts);
        Xs = fliplr(Xs);
        n = length(ts);
        Ls = zeros(3, 3, n);
        for i = 1:n
            Ls(:,:,i) = Ls_f(:,:,(n-i)+1);
        end
    else

        % Setup output arrays, n is the number of time steps to run
        % We include the first point in the output (hence the + 1)
        n = fix((tend - tstart)/dt) + 1;
        Xs = zeros(3, n);
        Ls = zeros(3, 3, n);
        ts = zeros(1, n);

        % Loop over time steps of path line
        Xi = X0;
        ti = tstart;
        % Store first point
        ts(1) = ti;
        Xs(:, 1) = Xi;
        Ls(:, :, 1) = AP_veloc_grad(Xi, ti, velocity_f, dx);
        for i = 2:n
            % Update point
            [Xi, ti] = integrator( Xi, ti, velocity_f, dt );
            % Store current point
            ts(i) = ti;
            Xs(:, i) = Xi;
            Ls(:, :, i) = AP_veloc_grad(Xi, ti, velocity_f, dx);
        end
    
    end

end

%
% Perform a single fourth-order Runge-Kutta step of timestep 
% dt along a pathline in an (possibly time varying) velocity
% field defined by the @velocity function starting from time
% t and location Xi. The position at the next step is given 
% by Xn
%
function  [Xn, tn] = AP_rk4_step( Xi, t, velocity_f, dt )

    k1 = velocity_f(Xi, t);
    k2 = velocity_f(Xi + 0.5*dt*k1 , t + 0.5*dt);
    k3 = velocity_f(Xi + 0.5*dt*k2 , t + 0.5*dt);
    k4 = velocity_f(Xi + dt*k3 , t + dt);
    
    tn = t + dt;
    Xn = Xi + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);

end

%
% Perform a single fourth-order Runge-Kutta step backwards in
% time
%
function  [Xn, tn] = AP_rk4_step_backwards( Xi, t, velocity_f, dt )

    k1 = -velocity_f(Xi, t);
    k2 = -velocity_f(Xi + 0.5*dt*k1 , t - 0.5*dt);
    k3 = -velocity_f(Xi + 0.5*dt*k2 , t - 0.5*dt);
    k4 = -velocity_f(Xi + dt*k3 , t - dt);
    
    tn = t - dt;
    Xn = Xi + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);

end

%
% Evaluate the velocity gradient tensor from a velocity field
% defined by the function velocity at point X and time t by
% centeral finite differencing with a fixed step size dx.
% I think this can be written in a cleaner way...
%
function [dVdX] = AP_veloc_grad( X, t, velocity_f, dx)

   dVdX = zeros(3,3);
   
   % Velocities at displaced positions
   Vpx = velocity_f(X+[ dx 0  0]', t);
   Vmx = velocity_f(X+[-dx 0  0]', t);
   Vpy = velocity_f(X+[0  dx  0]', t);
   Vmy = velocity_f(X+[0 -dx  0]', t);
   Vpz = velocity_f(X+[0  0  dx]', t);
   Vmz = velocity_f(X+[0  0 -dx]', t);
   
   % Derivatives
   dVdX(1, 1) = (Vpx(1) - Vmx(1)) / (2.0*dx);
   dVdX(2, 1) = (Vpx(2) - Vmx(2)) / (2.0*dx);
   dVdX(3, 1) = (Vpx(3) - Vmx(3)) / (2.0*dx);
   
   dVdX(1, 2) = (Vpy(1) - Vmy(1)) / (2.0*dx);
   dVdX(2, 2) = (Vpy(2) - Vmy(2)) / (2.0*dx);
   dVdX(3, 2) = (Vpy(3) - Vmy(3)) / (2.0*dx);

   dVdX(1, 3) = (Vpz(1) - Vmz(1)) / (2.0*dx);
   dVdX(2, 3) = (Vpz(2) - Vmz(2)) / (2.0*dx);
   dVdX(3, 3) = (Vpz(3) - Vmz(3)) / (2.0*dx);
end