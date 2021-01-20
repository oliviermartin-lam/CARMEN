classdef canaryTools
    
   methods (Static)
       
       
       function [alpha,delta] =  givesCoordinatesFromObject(object)

           object = upper(object);
           
           switch object
               
               case 'A47'
                   alpha  = [21,12,3.0];
                   delta = [38,36,45.0];
                                      
               case 'A453'
                   alpha  = [+1 36 44.400];
                   delta = [47 21 23.0];
                   
               case 'A396'
                   alpha  = [6 48 15.0];
                   delta = [41 05 44.0];
                   
               case 'A53'
                   alpha  = [+23 24 30.0];
                   delta = [40 53 55.0];
                   
               case 'A110'
                   alpha  = [+19 31 11];
                   delta = [34 56 02];
                   
               case 'A32'
                   alpha  = [+18 27 7.0];
                   delta = [26 52 34.0];
                   
               case 'A34'
                   alpha  = [+18 51 30.0];
                   delta = [10 20 3.0];
                   
               case 'A51'
                   alpha  = [+22 48 5];
                   delta = [39 17 3.0];
                   
               case 'A10'
                   alpha  = [+5 52 17.0];
                   delta = [32 34 36.0];
                   
               case 'AT1'
                   alpha  = [+13 41 38.2];
                   delta = [+7 36 21.0];
                   
               case 'ASTT1'
                   alpha  = [+13 41 38.2];
                   delta = [+7 36 21.0];
                   
               case 'A12'
                   alpha  = [+6 1 9.0];
                   delta = [23 20 29.0];
               
               otherwise  
                   delta = [];
                   alpha = [];
           end
       end


       function out = airmassFromDate(date,alphaJ2000,deltaJ2000)
           
           %extract date from path
           Y = str2double(date(1:4));
           M = str2double(date(5:6));
           D = str2double(date(7:8));
           H = str2double(date(10:11))-1.;%correction temps local vers temps UTC
           Mm = str2double(date(13:14));
           S = str2double(date(16:17));
           %defines coordinates
           alpha = canaryTools.alpha2deg(alphaJ2000(1),alphaJ2000(2),alphaJ2000(3));
           delta = canaryTools.dec2deg(deltaJ2000(1),deltaJ2000(2),deltaJ2000(3));
           
           out = canaryTools.airmassFromCoordinates(Y, M, D, H, Mm, S, alpha, delta);
       end
       
       function out =  alpha2deg(h, m, s)
           out = canaryTools.dec2deg(h, m, s)*15;
       end
       
       function out =  dec2deg(h, m, s)
           out = h+m/60.+s/3600.;
       end
       
       function out =  airmassFromCoordinates(Y, M, D, H, Mm, S, alpha, delta, varargin)
           inputs = inputParser;
           inputs.addRequired('Y',@isnumeric);           
           inputs.addRequired('M',@isnumeric);       
           inputs.addRequired('D',@isnumeric);           
           inputs.addRequired('H',@isnumeric);  
           inputs.addRequired('Mm',@isnumeric);           
           inputs.addRequired('S',@isnumeric);  
           inputs.addRequired('alpha',@isnumeric);           
           inputs.addRequired('delta',@isnumeric);  
           inputs.addParameter('latitude',28.7567,@isnumeric);
           inputs.addParameter('longitude',-17.8917,@islogical);
           inputs.parse(Y, M, D, H, Mm, S, alpha, delta,varargin{:});
           
           latitude  = inputs.Results.latitude;
           longitude = inputs.Results.longitude;
  
           h = (H-12)/24. + Mm/24./60 + S/24./3600;
           J = (1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075;
           JJ= double(J) + h;            
           D = JJ - 2451545.;
           
           GMST = 18.697374558 + 24.06570982441908 * D;
           lst = GMST+longitude/15.;
           lst = mod(lst,24);
            
           H = lst - alpha;
           r = pi/180;
           sina = sin(latitude*r)*sin(delta*r)+cos(latitude*r)*cos(delta*r)*cos(H*r);                      
           zenithAngle = pi/2.-asin(sina);
           out =  1/cos(zenithAngle);
       end                                 

        %%
       function out =  getValidSubapArray( nLenslet, rext, rint )
           % DOCUMENT valid = getValidSubapArray( nLenslet, rext, rint )
           
           %Example :
           %> getValidSubapArray( 5, 1, 0.285 )
           %[0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0]
           
   
 
           
           % tip-tilt sensor only.
           if nLenslet==1
               out = [1];
           else
               % to avoid some bug that eliminates useful central subapertures when obs=0.286
               if nLenslet==7 && (rint>0.285 && rint<0.29)
                   rint=0.285;
               end
                   
               x     = linspace(-1+1./nLenslet,1-1./nLenslet,nLenslet);               
               [x,y] = meshgrid(x);
               r     = hypot(x,y);
               out   = (r<=rext & r>=rint);
               out   = out(:);               
           end
       end


       function out = mirror_SH7(s_orig, sym, varargin)
           % DOCUMENT smir = mirror_SH7(s_orig, sym, t=, inverse=)      
           inputs = inputParser;
           inputs.addRequired('s_orig',@isnumeric);           
           inputs.addRequired('sym',@isnumeric);       
           inputs.addParameter('flagTranspose',0,@isnumeric);
           inputs.addParameter('flagInverse',false,@islogical);
           inputs.parse(s_orig,sym,varargin{:});
           
           flagTranspose = inputs.Results.flagTranspose;
           flagInverse = inputs.Results.flagInverse;
           out = canaryTools.mirror_SH_data(s_orig,7,1.0,0.285,sym,'flagTranspose',flagTranspose,'flagInverse',flagInverse);
       end
       
       function out = mirror_SH_data(s_orig, nLenslets, radius, telObs, sym, varargin)
           % DOCUMENT smir = mirror_SH_data(s_orig, nssp, ext_rad, obst_rad, sym, t=, inverse=)
           
           %<s_orig> is an array of slopes, or an interaction matrix
           %<nssp> is the number of subapertures (7)
           %<ext_rad> is the pupil external radius (1.07)
           %<obst_rad> is the radius of central obstruction
           %<sym> is the symmetry number, between 0 to 7, i.e.
           %      0 = nothing
           %      1 = symmetry around a vertical axis
           %      2 = symm around an horiz axis
           %      4 = symm around axis x=y (transpose)
           %For several symmetries at once, just add.
           
           %<t> allows to work on data that were transposed, i.e. where the dim
           %    of #subap is placed in 2nd position, like in interaction matrices.
           %<inverse>=1 will apply the reciprocal transformation
           
           inputs = inputParser;
           inputs.addRequired('s_orig',@isnumeric);       
           inputs.addRequired('nLenslets',@isnumeric );                        
           inputs.addRequired('radius',@isnumeric);       
           inputs.addRequired('telObs',@isnumeric);       
           inputs.addRequired('sym',@isnumeric);       
           inputs.addParameter('flagTranspose',0,@isnumeric);
           inputs.addParameter('flagInverse',[],@islogical);
           inputs.parse(s_orig,nLenslets,radius,telObs,sym,varargin{:});
           
           flagTranspose = inputs.Results.flagTranspose;
           flagInverse = inputs.Results.flagInverse;

           
           
           sym = mod(sym,8);
           
           if flagTranspose
               if sym==5
                   sym=6;
               elseif sym==6
                   sym=5;
               end
           end
                                          
           if flagTranspose==0
               out = s_orig;
           else
               s       = transpose(s_orig);
               valid   = canaryTools.getValidSubapArray( nLenslets, radius, telObs );
               nn      = find(valid);
               msk     = zeros(nLenslets);
               msk(nn) = 1:length(nn);
               sm      = s;
               xs      = 1;
               ys      = 1;
               
               if sym ==1
                   msk = msk(1:end-1,:);
                   xs=-1;
               end
               
               if sym ==2
                   msk = msk(:,1:end-1);
                   ys  = ys-1;
               end
               if sym ==4
                   msk              = msk';
                   ind              = msk(nn);
                   nssp             = length(ind);
                   sm(1:nssp,:)     = s(ind+nssp,:) * ys;
                   sm(nssp+1:end,:) = s(ind,:) * xs;
               else
                   ind              = msk(nn);
                   nssp             = numberof(ind);
                   sm(1:nssp,:)     = s(ind,:) * xs;
                   sm(nssp+1:end,:) = s(ind+nssp,:) * ys;
               end
               
               if flagTranspose==0
                   out= sm;
               else
                   out= sm';
               end
           end
       end              
   end
end

