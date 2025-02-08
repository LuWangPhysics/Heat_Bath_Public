        function obj=init_field(obj,C,E_iter)
      
        obj.E0 = E_iter;                                                    % peak electric field strength     
        obj.lambda_si=3200e-9;   
        obj.om = C.c*(2*pi/obj.lambda_si);                                 % circular frequency; 0.057 [at.u.] approx 800nm


        obj.tau=20e-15;                                                   %pulse duration 
        obj.type_N=16000;
      
              % N_mesh=8;
              % t_lim=11*obj.tau;


        % %%for loose mesh
             t_lim=6.5*obj.tau;
              N_mesh=1.1;
 
  
        obj.t = linspace(-t_lim,t_lim,N_mesh*8e3+1);                       % time vector                                          
        obj.dt = obj.t(2)-obj.t(1);    
        obj.t_win_effect=fix(t_lim/2.1/obj.dt);                          %effetive time windwo for ionizaiton 
        
        obj.dw=2*pi/obj.dt/length(obj.t);
        w_limit=length(obj.t).*obj.dw/2;
        obj.w=linspace(-w_limit,w_limit,length(obj.t));
        
        obj.Et = obj.E0*exp(-(obj.t/obj.tau).^2).*cos(obj.om*obj.t);       % field function of time
        obj.At = -cumtrapz(obj.Et)*obj.dt;                                   % vector potential 


        obj.t_start=fix(length(obj.t)/2)-obj.t_win_effect;
        obj.t_end=fix(length(obj.t)/2)+obj.t_win_effect;
        end 

