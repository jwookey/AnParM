% AP_update_finite_strain - update the finite strain tensor based on an applied 
%                           velocity gradient..
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage: [FST_new,c,R,ph] = AP_update_finite_strain(FST,vgrad,t,PrevR)
%
%     Inputs: 
%        FST    : Previous finite strain tensor (3x3)
%        vgrad  : velocity gradient tensor (3x3)
%        t      : deformation time
%        PrevR  : the previous rotation matrix (3x3)*
%    
%     Outputs:
%        FST_new : updated finite-strain tensor
%        c       : axis lengths of the (new) finite strain ellipse
%        R       : rotation matrix required to map a vector/tensor in
%                  the original cartesian reference frame onto that
%                  of the deformed FSE.
%        ph      : minimum angle of this rotation
%    
%      *The code calculates the principle axes of the finite strain ellipse, and 
%      calculates the arrangement of these axes which imply the smallest rotation
%      from the previous rotation. This is potentially rather a rate-limiting step.
%
% Reference.
% ~~~~~~~~~~
%
%  McKenzie 1979
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


function [FST_new,c,R,ph] = AP_update_finite_strain(FST,vgrad,t,PrevR) ;
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