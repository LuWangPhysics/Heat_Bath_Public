classdef Rho
properties 
    rho_0;
    rho_1;
    rho_2;

    J_t;
    int_t1_arr;
    hbar;
    q;
    rho_2_no_heat;
    J_E_t_tau;


    N_rho;
    N_t;
    Norm_k;
    rho_nm;
    rho_k;
    rho_nn;
   magic_operation;

end
  methods
          function obj=init(obj,J_t,MF,C,k_select)
              obj.Norm_k=length(k_select);
              obj.N_rho= MF.t_end-MF.t_start+1;
              if mod(length(MF.t),2)==0
                  obj.N_t=fix(length(MF.t)/2)+1; 
              else
                  obj.N_t=fix(length(MF.t)/2)+2;
              end
              obj.rho_0=zeros(1,obj.N_rho);
              obj.rho_1=zeros(1,obj.N_rho);
              obj.rho_2=zeros(1,obj.N_rho);
              obj.rho_2_no_heat=zeros(1,obj.N_rho);
              obj.rho_nm=0;
              obj.rho_nn=0;
              obj.J_t=J_t;
 
              obj.hbar=C.hbar;
              obj.q=C.q;
              %the first dimension is the one contracted
        
              if length(MF.t)<MF.type_N
                      obj.J_E_t_tau=diag(ones(1,obj.N_rho).*0.5);
                       % -------------------------------------------------------
                       %construct 2D matrix for nested t1, t2 integral with t1>t2 
                       %-J_E_t_tau(:,1) is t1-t2, the array need to be fliped when multiply with the f_1(t2)
                       %-J_E_t_tau(1,:) is t1
                       % -------------------------------------------------------
                        for iter_tau=MF.t_start+1:MF.t_end
                          temp=obj.J_t(obj.N_t:(obj.N_t+iter_tau-MF.t_start-1));
                          obj.J_E_t_tau(1:(iter_tau-MF.t_start),iter_tau-MF.t_start+1)=flip(temp);

                        end
                        obj.magic_operation=@(f1_HI,MF) obj.manipulate_short(f1_HI,MF);
              else
               %for fine dt data points, to save memory, it goes back to
               %brute force integration
                        obj.magic_operation=@(f1_HI,MF) obj.manipulate_long(f1_HI,MF);
              end
   
           
          end

          function output=manipulate_short(obj,f1_HI,MF)
             output =(f1_HI*(obj.J_E_t_tau)).*MF.dt;
          end


          function output=manipulate_long(obj,f1_HI,MF)

                %calculate matrix multiplication for each sub sections
                    output=zeros(size(f1_HI));
                     output(1)=f1_HI(1)*1/2;
                        for iter_tau=MF.t_start+1:MF.t_end-1
                          temp1=obj.J_t(obj.N_t:(obj.N_t+iter_tau-MF.t_start-1));
                          temp2=flip([1/2,temp1]);
                          output(iter_tau-MF.t_start+1)=sum(f1_HI(1:length(temp2)).*temp2).*MF.dt;
                        end


          end


          function obj=rho_update(obj,MF,C,jo)
             t_range=MF.t_start:MF.t_end;
             V1=sqrt((MF.Eps_total(t_range)+MF.Eps_0(t_range))./(2*MF.Eps_total(t_range)));
             %---------------------------------------
             %define the f_1 functin for integration
             %---------------------------------------
             temp1=(MF.Omega_dipole.*MF.Eps_0_dot-MF.Eps_0.*MF.Omega_dipole_dot)...
                 .*(MF.Eps_total-MF.Eps_0)./(2.*MF.Eps_total.^3);
             temp2=1i.*MF.Omega_dipole./C.hbar./2;
         
             f1=temp1(t_range)-temp2(t_range);

             %---------------------------------------
             %perform the integration
             %---------------------------------------

             %make sure the right dimension is multiplied
             int_eps_total=cumtrapz(MF.Eps_total(t_range));
             f1_HI=f1.*exp(-1i*(MF.dt*int_eps_total./C.hbar));
             first_int=obj.magic_operation(f1_HI,MF);

             obj.rho_k=2.*V1.^2.*real(cumtrapz(conj(f1_HI).*first_int.*MF.dt))./obj.Norm_k;
             obj.rho_2=obj.rho_2+obj.rho_k;

             %---------------------------------------------
             %no heat bath
             %---------------------------------------------          


             obj.rho_2_no_heat=obj.rho_2_no_heat+V1.^2.*abs(cumtrapz(f1_HI).*MF.dt).^2./obj.Norm_k;
            

          end


    end
end