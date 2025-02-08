     function obj=init_material(obj,m_name,C)

          switch m_name
              case 'ZnO'
                obj.eg = 0.129.*C.Energy_si;                                       % min bandgap for ZnO
                obj.delx = 0.17*C.Energy_si;                                                   % bandwidth along direction x, ZnO

                obj.ax = 5.32.*C.r_si;                                             % length of unit cell along direction x, ZnO
                obj.d0 = 3.64*C.r_si;                                                     % dipole d(kx=0) at gamma point, ZnO
              case 'SiO2'
    
                % ---------------------------------------------------------- model parameters dielectrics --------------------------------------------
                obj.eg = 0.32.*C.Energy_si;
                obj.delx =[0.065,-0.0035,-0.001,-0.0007] *C.Energy_si;
  
                obj.ax = 9*C.r_si;
                obj.d0 = 6.5*C.r_si;    
                % dipole along kx in K.P approximation ==> d(0) * (Eb(0)/Eb(kx))
              otherwise
                   disp('Human! You have not defined me yet!')
    
          end


                                                   
           obj.k = linspace(-pi/obj.ax,pi/obj.ax,101);                     % crystal momentum vector along direction x                                                                             
 
            
  
        end
