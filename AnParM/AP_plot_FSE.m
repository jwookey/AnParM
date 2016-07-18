function [] = AP_plot_FSE(FSE_axes,RM,varargin)

%  ** Set defaults, these can be overriden in the function call
      i3d = 0 ;
      planeSel = 'max' ;
      
%  ** process the optional arguments
      iarg = 1 ;
      while iarg <= (length(varargin))
         switch lower(varargin{iarg})
            case '3d'
               i3d = 1 ;
               iarg = iarg + 1 ;
            case 'plane'
               planeSel = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            case 'disp'
               dispvect = varargin{iarg+1} ;
               iarg = iarg + 2 ;
            otherwise 
               error(['Unknown option: ' varargin{iarg}]) ;   
         end   
      end

      % if i3d & ~isNaN(planeSel)
      %    error('3d plotting and plane select are incompatble options')
      % end   

      if ischar(planeSel) 
         % check that it is acceptable
         switch lower(planeSel)
         case '12'
         case '13'
         case '23'
         case 'max'   
         otherwise 
            error(['Invalid plane selection: ' planeSel]) ;   
         end
      else
         error('Arbitrary planes are currently not supported.')
      end
         
%  ** for 2D plotting
   if i3d
      error('3D plotting not supported.') ;
   else
      AP_plot_FSE_2D(FSE_axes,RM,planeSel)
   end

   daspect([1 1 1]) ;

end

function [] = AP_plot_FSE_2D(FSE_axes,RM,planeSel)
   
   % rotate axis vector
   %RFSE = RM*FSE_axes

   XYZ=[cos([-pi:pi/101:pi]).*FSE_axes(1) ; ...
        sin([-pi:pi/101:pi]).*FSE_axes(2) ; ...
        zeros(1,length([-pi:pi/101:pi])) ] ;
   
   
   %figure
   
   %plot(XYZ(1,:),XYZ(2,:),'k-') ;
   
   XYZR = RM*XYZ ;
   hold on ;
   plot(XYZR(1,:),XYZR(2,:),'k-') ; 
   
   

end
