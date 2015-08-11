% Test harness to begin conversion of odfcalc   
function AP_odfcalc_Ol_simple()
      
   % hardwire inputs
   
   % velocity gradient tensor
   vgrad = [-1 0 0 ; 0 0.5 0 ; 0 0 0.5] ;
   
   % power law exponent
   rn = 3.5 ;
   % CRSS
   tau = [0.3333 0.6667 1.0] ;
   
   % number of grains
   ngrains = 2000 ;
   
   % t
   dt = 0.6 ;
   
   % only one step allowed currently
   
   % check vgrad trace
   
   %
   p12 = log(tau(1)/tau(2)) ;
   p23 = log(tau(2)/tau(3)) ; 
   p13 = log(tau(1)/tau(3)) ;
   
   nslip = 3 ;
   
   eps = zeros([3 3 3]) ;
   
   eps(1,2,3) = 1 ;
   eps(2,3,1) = 1 ;
   eps(3,1,2) = 1 ;
   eps(2,1,3) = -1 ;
   eps(3,2,1) = -1 ;
   eps(1,3,2) = -1 ;
   
   % build random texture
   [ eulers ] = MVT_make_random_texture( ngrains ) ;
   
   % make it an n x 3 list, and convert to radians
   
   eulers = (pi/180) * eulers' ;
      
   % start of calculations
   
   %eulers = load('my_eulers.dat')*pi/180 ;
   
   % Strain Rate Tensor
   bige = zeros(3,3) ;
   for i = 1:3
      for j = 1:3
         bige(i,j) = (vgrad(i,j) + vgrad(j,i))/2.d0 ;
      end
   end
   
   % rotation rates vector
   bigom = [(vgrad(3,2) - vgrad(2,3))/2 ...
            (vgrad(1,3) - vgrad(3,1))/2 ...
            (vgrad(2,1) - vgrad(1,2))/2] ;
   
   
   % initial time
   t = 0 ;
   
   % initial FSE parameters
   r12 = 0 ; r23 = 0 ; r13 = 0 ; % spherical
   
   [bigA1,bigA2,bigA3,bigB1,bigB2,bigB3] = ...
      amplitudes(p12,p23,p13,r12,r23,r13) ;

   % loop over grains
   for ig=1:ngrains
      % calculate the spin
      a=dircos(eulers(ig,:)) ;
      [en,el] = slipvectors(a) ;
      
      %
      gamdot = sliprates(nslip,en,el, bige, ...
                         bigA1,bigA2,bigA3, ...
                         bigB1,bigB2,bigB3) ;
      omc = crystalspin(eps,en,el,gamdot,nslip) ;
      
      % total spin
      om = omc + bigom ;
      
      % spin rate
      gdot = gdotcalc(om,eulers(ig,:)) ;
      
      
      % apply the spin rate for the appropriate time
      eulers(ig,:) = eulers(ig,:) + dt.*gdot ;
      disp( eulers(ig,:) ) ;
      
   end      
   


   eulers_out = 180/pi * eulers' ;
   MVT_write_VPSC_file('simple.out', eulers_out, 'Test output')
   
         
end

function [gdot] = gdotcalc(om,eul)
   cosph = cos(eul(1)) ;
   sinph = sin(eul(1)) ;
   costh = cos(eul(2)) ;  
   sinth = sin(eul(2)) ;         
   cotth = costh/sinth ;
   
   % phdot, thdot, psdot
   gdot = [om(3) + cotth*(om(2)*cosph - om(1)*sinph) ...
           om(1)*cosph + om(2)*sinph ...
           (om(1)*sinph - om(2)*cosph)/sinth ] ;
   
end

function [omc] = crystalspin(eps,en,el,gamdot,nslip)
   omc = [0 0 0] ;
   for i=1:3
      for is=1:nslip
         for j=1:3
            for k=1:3
               omc(i) = omc(i) - 0.5*eps(i,j,k)*en(is,j)*el(is,k)*gamdot(is) ;
            end
         end
      end
   end
end
               

function [gamdot] = sliprates(nslip,en,el, bige, ...
                              bigA1,bigA2,bigA3, ...
                              bigB1,bigB2,bigB3) ;
   bigS=zeros(3,3) ;
   gamdot = [0 0 0] ;
   for is=1:nslip
      % calculate the Schmidt tensor
      for i=1:3
         for j=1:3
            bigS(j,i) = 0.5*(en(is,i)*el(is,j) + en(is,j)*el(is,i)) ;
         end
      end
      gamdot(is) = ...
         bigA1(is)*bige(1,1)*bigS(1,1) ...
           + bigA2(is)*(bige(2,2)*bigS(1,1) ...
           + bige(1,1)*bigS(2,2)) ...
           + bigA3(is)*bige(2,2)*bigS(2,2) ...
           + 4.0*bigB1(is)*bige(1,2)*bigS(1,2) ...
           + 4.0*bigB2(is)*bige(2,3)*bigS(2,3) ...
           + 4.0*bigB3(is)*bige(1,3)*bigS(1,3) ;
   end
   
   
   
end

function [a] = dircos(eul) ;
   a=zeros(3,3) ;
   a(1,1) = cos(eul(1))*cos(eul(3)) - cos(eul(2))*sin(eul(1))*sin(eul(3)) ;
   a(1,2) = cos(eul(3))*sin(eul(1)) + cos(eul(1))*cos(eul(2))*sin(eul(3)) ;
   a(1,3) = sin(eul(3))*sin(eul(2)) ;
   a(2,1) = - cos(eul(3))*cos(eul(2))*sin(eul(1)) - cos(eul(1))*sin(eul(3)) ;
   a(2,2) = cos(eul(1))*cos(eul(3))*cos(eul(2)) - sin(eul(1))*sin(eul(3)) ;
   a(2,3) = cos(eul(3))*sin(eul(2)) ;
   a(3,1) = sin(eul(1))*sin(eul(2)) ;
   a(3,2) = - cos(eul(1))*sin(eul(2)) ;
   a(3,3) = cos(eul(2)) ;
end

function [en,el] = slipvectors(a)
   en = zeros(3,3) ;
   el = zeros(3,3) ;
   
   for i=1:3
      %  Slip system (010)[100]
      en(1,i) = a(2,i) ;
      el(1,i) = a(1,i) ;
      %  Slip system (001)[100]
      en(2,i) = a(3,i) ;
      el(2,i) = a(1,i) ;
      %  Slip system (010)[001]
      en(3,i) = a(2,i) ;
      el(3,i) = a(3,i) ;
   end
end
   
function [bigA1,bigA2,bigA3,bigB1,bigB2,bigB3] = ...
   amplitudes(p12,p23,p13,r12,r23,r13)
   
   r0 = (2./3.)*sqrt(r12^2 + r12*r23 + r23^2) ;
   
   calA11 = 2.241 + 0.3993*r23^2 + 1.104*r12*r23 ...
          + 1.104*r12^2 + 0.7619*r23^4 + 2.507*(r23^3)*r12 ...
          + 5.518*(r23^2)*r12^2 + 6.023*r23*r12^3 + 3.012*r12^4 ;
   
   calA12	 = 2.241 + 0.3993*r12^2 + 1.104*r12*r23 ...
          + 1.104*r23^2 + 0.7619*r12^4 + 2.507*(r12^3)*r23 ...
          + 5.518*(r23^2)*r12^2 + 6.023*r12*r23^3 + 3.012*r23^4 ;
   
   calA13 = 2.241 + 0.3993*r13^2 - 1.104*r12*r13 ...
          + 1.104*r12^2 + 0.7619*r13^4 - 2.507*(r13^3)*r12 ...
          + 5.518*(r13^2)*r12^2 - 6.023*r13*r12^3 + 3.012*r12^4 ;
   
   calB11 = 1.662 + 0.2046*r23^2 + 0.1992*r12*r23 ...
          - 0.7517*r12^2 - 0.01853*r23^4 - 0.02831*(r23^3)*r12 ...
          - 0.4396*(r23^2)*r12^2 - 0.4246*r23*r12^3 ...
          + 0.2085*r12^4 ;
   
   calB12 = 1.662 + 0.2046*r12^2 + 0.1992*r12*r23 ...
          - 0.7517*r23^2 - 0.01853*r12^4 - 0.02831*(r12^3)*r23 ...
          - 0.4396*(r23^2)*r12^2 - 0.4246*r12*r23^3 ...
          + 0.2085*r23^4 ;
   
   calB13 = 1.662 + 0.2046*r23^2 - 0.1992*r13*r23 ...
          - 0.7517*r13^2 - 0.01853*r23^4 + 0.02831*(r23^3)*r13 ...
          - 0.4396*(r23^2)*r13^2 + 0.4246*r23*r13^3 ...
          + 0.2085*r13^4 ;
   
   cala1 = 1.000 - 0.0295*p12 - 0.0130*p23 ...
         - 0.00347*p12*p23 - 0.00743*p12^2 - 0.00333*p23^2 ;
   cala2 = 1.000 + 0.0295*p12 - 0.0130*p13 ...
         + 0.00347*p12*p13 - 0.00743*p12^2 - 0.00333*p13^2 ;
   cala3 = 1.000 + 0.0295*p23 + 0.0130*p12 ...
         - 0.00347*p12*p23 - 0.00743*p23^2 - 0.00333*p12^2 ;
   
   bigA1 = [cala1*(2.0*calA11 + 2.0*calA12 - calA13) ...
            cala2*(2.0*calA11 + 2.0*calA12 - calA13) ...
            cala3*(2.0*calA11 + 2.0*calA12 - calA13) ] ;
   
   bigA2 = [-0.5*cala1*(calA11 - 5.0*calA12 + calA13) ...
            -0.5*cala2*(calA11 - 5.0*calA12 + calA13) ...
            -0.5*cala3*(calA11 - 5.0*calA12 + calA13) ] ;
   
   bigA3 = [-1.0*cala1*(calA11 - 2.0*calA12 - 2.0*calA13) ...
            -1.0*cala2*(calA11 - 2.0*calA12 - 2.0*calA13) ...
            -1.0*cala3*(calA11 - 2.0*calA12 - 2.0*calA13) ] ;
   
   bigB1 = [cala1*calB11 cala2*calB11 cala3*calB11] ;
   
   bigB2 = [cala1*calB12 cala2*calB12 cala3*calB12] ;
   
   bigB3 = [cala1*calB13 cala2*calB13 cala3*calB13] ;

   end
   
   