classdef Material_and_Field
    properties
        %--------------------
        %field parameters
        %--------------------
        dt;
        t;
        om;
        dw;
        w;
        tau;
        E0;
        lambda_si;
        t_win_effect;
        t_start;
        t_end;

        Et;
        At;
        type_N;
        %--------------------
        %material
        %--------------------
        k;
        eg;
        delx;
        ax;
        d0;
        Eps_total;
        Omega_dipole;
        Omega_dipole_dot;
        Eps_0;
        Eps_0_no_t;
        Eps_0_dot;
        vx;
        d;

    end

    methods


   
        function obj=renew_kt(obj,C,k_iter)

       Kt = obj.k(k_iter)+C.e.*obj.At./C.hbar;                                      %the k after shifted by A                          % kt = K-q A(t)/hbar 
 
       %construct eps0 and v
        obj.Eps_0=obj.eg;
        obj.Eps_0_no_t=obj.eg;
        obj.vx=0;
        for i_iter=1:length(obj.delx)
                 obj.Eps_0=obj.Eps_0+obj.delx(i_iter).*(1-cos(Kt.*i_iter.*obj.ax));
                 obj.Eps_0_no_t=obj.Eps_0_no_t+obj.delx(i_iter).*(1-cos(obj.k.*i_iter.*obj.ax));
                 obj.vx=obj.vx+obj.ax.*i_iter.*obj.delx(i_iter).*sin(Kt.*i_iter.*obj.ax)./C.hbar;
        end

        obj.d=obj.eg.*obj.d0./obj.Eps_0;
        obj.Omega_dipole=2*C.q.*obj.Et.*obj.d;                               
        obj.Eps_total=  sqrt( obj.Eps_0.^2 +obj.Omega_dipole.^2 );
        obj.Eps_0_dot=obj.vx.*C.q.*obj.Et;
        obj.Omega_dipole_dot=diff([obj.Omega_dipole,obj.Omega_dipole(end)])./obj.dt;
        end

    end

end 