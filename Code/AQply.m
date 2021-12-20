classdef AQply < handle
   properties
       %% Base Material Properties (Undamaged)
       E1
       E2
       G12
       nu12
       nu21
%        alph1
%        alph2
       %% Damage
       df  % fiber damage
       dm  % matrix damage
       ds  % shear damage
       D   % damage coeff
       %% Material Properties (Damage Applied)
%        d_E1
%        d_E2
%        d_G12
%        d_nu12
%        d_nu21
       %% Stress-Strain Relation Matricies
       Q126   % Ply coordinates
       Qxys   % Laminate coordinates
       %% Laminate Positioning
       t   % thickness
       ori % theta orientation angle
       z1  % bottom stacking position
       z2  % top stacking position
       T   % transformation matrix
       %% Results
       epsxys  % laminate strains
       eps126  % ply orientation strains
       sig126  % ply stresses
   end
   
   methods

       %initialize
       function obj = AQply(E1,E2,G12,nu12,nu21)
           % send a vector of material properties
           obj.E1 = E1;
           obj.E2 = E2;
           obj.G12 = G12;
           obj.nu12 = nu12;
           obj.nu21 = nu21;
%            obj.alph1 = alph1;
%            obj.alph2 = alph2;
       end
                 
       function calcDmgProps(obj, df, dm)
           %Outlined in Lapcyzk, Hurtado (2007)
           obj.df = df;
           obj.dm = dm;
           dft = df;  % tensile fiber damage
           dfc = df;  % compressive fiber damage
           dmt = dm;  % tensile matrix damage
           dmc = 0;   % compressive matrix damage
           obj.ds = 1-(1-dft)*(1-dfc)*(1-dmt)*(1-dmc);
           obj.D = 1-(1-df)*(1-dm)*obj.nu12*obj.nu21;
%            obj.d_E1 = obj.E1 *(1-obj.d(1));
%            obj.d_E2 = obj.E2 *(1-obj.d(2));
%            obj.d_G12 = obj.G12 *(1-obj.d(3));
%            obj.d_nu12 = obj.nu12 *(1-obj.d(4));
%            obj.d_nu21 = obj.nu21 *(1-obj.d(5));
       end
       
       function makePlyQ(obj)
%            Q11 = obj.d_E1 / (1-obj.d_nu12*obj.d_nu21);
%            Q22 = obj.d_E2 / (1-obj.d_nu12*obj.d_nu21);
%            Q12 = (obj.d_nu21*obj.d_E1) / (1-obj.d_nu12*obj.d_nu21);
%            Q66 = obj.d_G12;
           
           Q11 = 1/obj.D * (1-obj.df) * obj.E1;
           Q22 = 1/obj.D * (1-obj.dm) * obj.E2;
           Q12 = 1/obj.D * (1-obj.df)*(1-obj.dm)*obj.nu21*obj.E1;
           Q66 = (1-obj.ds)*obj.G12;
           
           obj.Q126 = [...
               Q11, Q12,   0 ; ...
               Q12, Q22,   0 ; ...
                0 ,  0 , Q66];
       end
       
       function placePly(obj,thetaDeg,thick,zbot)
           obj.ori = thetaDeg;
           obj.t = thick;
           obj.z1 = zbot;
           obj.z2 = zbot+thick;
       end
       
       function offsetZvals(obj, offset)
           %offset should be 1/2 the laminate thickness
           obj.z1 = obj.z1 - offset;
           obj.z2 = obj.z2 - offset;
       end
       
       function makeTmatrix(obj)
           m = cos(pi/180*obj.ori);
           n = sin(pi/180*obj.ori);
           
           obj.T = [...
               m^2 , n^2 , 2*m*n ; ...
               n^2 , m^2 ,-2*m*n ; ...
              -m*n , m*n , m^2 - n^2 ];
       end
        
       function calcQxys(obj)
           % See CLT.  Some shear terms have coeffiecents
           modQ126 = obj.Q126.*[1 1 1; 1 1 1; 1 1 2] ;
           modQxys = obj.T \ modQ126 * obj.T;
           obj.Qxys = modQxys * diag([1,1,0.5]);
       end
       
       function applyxysStrain(obj, xysStrainVec)
           %% epsx, epsy, gammas
           %% eps1, eps2, ½*gamma6
           obj.epsxys = xysStrainVec(1:3);
           if any( xysStrainVec(4:6) > 1e-10 )
               sprintf('WARNING: Non-zero Curvatures!!!')
           end
           obj.eps126 = diag([1,1,2])*(obj.T*(diag([1,1,0.5])*obj.epsxys));  
           obj.sig126 = obj.Q126 * obj.eps126;          
       end
       
       function Qat = calc_alphabar(obj, barangle)
           Angle between ply fibers and the "90 deg" ply according to
           Varna et al         
           theta = barangle - obj.ori;  
           m = cos(pi/180*theta);
           n = sin(pi/180*theta);
           
           Tk = [...
               m^2 , n^2 , 2*m*n ; ...
               n^2 , m^2 ,-2*m*n ; ...
              -m*n , m*n , m^2 - n^2 ];
          
           abar = Tk' * [obj.alph1 ; obj.alph2 ; 0];       
           Qkbar = Tk' * obj.Q126;           
           Qat = Qkbar * abar *obj.t;
       end
   end
end
           
           
           
       
           
           