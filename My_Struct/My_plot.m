classdef My_plot
    properties
    end
    methods
        function plot_E_field(obj,E)


            figure('Name','laser electric field','NumberTitle','off');
            plot(E.t./1e-15,E.Et)
            xlabel('Time (fs)')
            ylabel('E field strengh (V/m)')
        end

        function plot_rho_2(obj,E,R,T)

            figure
            plot(E.t(E.t_start:(E.t_start+length(R.rho_2)-1))./1e-15,[R.rho_2_no_heat-R.rho_2_no_heat(1);R.rho_2-R.rho_2(1)])
            set(gca, 'YScale', 'log')
            xlim([-4*E.tau*1e15,4*E.tau*1e15])
            legend('no heat bath','heat bath')
            xlabel('Time (fs)')
            ylabel('Ionization rate')
            title(['Temperature=' num2str(T) 'K'])
        end

        function J_t_int_FT_compare(obj,Jt_FT,Jt_int,MF,name_str,C,j0)
             name_select=[name_str,'_T=',num2str(C.T),'_Jo=' num2str(j0)];
            N_half=fix(length(MF.t)/2);
            plot_range =(N_half):(N_half+MF.t_win_effect);    
            figure('Name',name_select,'NumberTitle','off');
            subplot(1,2,1)
            plot(MF.t(plot_range)./1e-15,[real(Jt_FT(plot_range));imag(Jt_FT(plot_range))],'LineWidth',1,'Marker','square')
            xlim([0,5])
            xlabel('fs')
            ylim([-1e-3,1])
            title('numerical FFT')

            subplot(1,2,2)
            plot(MF.t(plot_range)./1e-15,[real(Jt_int(plot_range));imag(Jt_int(plot_range))],'LineWidth',1,'Marker','square')
            xlim([0,5])
            ylim([-1e-3,1])
            xlabel('fs')
            title('numerical integration')
    

        end
    end 
end