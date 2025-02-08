%Lulu's best code ever!
clear all;
restoredefaultpath;
addpath('My_Struct');
addpath('My_Functions');
addpath('My_Output');

% ========================== define constants==============================
P=My_plot;
C=Const;
C=C.init();
% ========================== material and driving laser parameters ========
MF=Material_and_Field;
M_type='ZnO';
MF=MF.init_material(M_type,C);                                              % choose between 'Dielectric','ZnO'

                            


% ========================== set-up heat bath general T numerical =========
C.T=300;                                                                % temperature K
jo=5;
E_in=1.5e9;
om_cutoff =200*(2*pi*1e12);                                                  % cutoff frequency
       
MF=MF.init_field(C,E_in);

        
% ========================== set-up heat bath general T numerical =========
%choose name among 'Dybe','OM','Gauss' 'Shift_Gauss','under_damp'
name_select='OM';
[Jt_FT,Jt_int]=numerical_general_heat_bath(C,MF,om_cutoff,jo,name_select);

% =================== set-up matrix for convolution integral for rho_2_22 =


k_select=1:length(MF.k);

R=Rho;
R=R.init(Jt_FT,MF,C,k_select);


for k_iter=k_select
        %----------------------------------------
        %renew the time dependent material 
        %----------------------------------------
        MF=MF.renew_kt(C,k_iter);
        
        %----------------------------------------
        %calculate the ionization for a given k
        %----------------------------------------
        R=R.rho_update(MF,C,jo);
end


P.plot_rho_2(MF,R,C.T);

