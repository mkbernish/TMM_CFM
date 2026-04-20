function gQpm1 = allocP(P,Qp)
Pmax=0.00320513285659728;
OT=0.00863364097132997;
m=3.79146798299876E-19   ;      %(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
I = 100;
Pchl=Pmax*(1-exp(-OT*I)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
Bchl = m/Pchl;
E=0.7742;
Achl = (1+E)/Pchl;
%========================
% Growth rate under P
%========================
Ap_RNA = 4.23e-3;Apho = 1.6e1;
Abio = 2.71e-1;Y_np_RNA = 3.8/1;
Qc_Pro_Other = 2.4e-1;
Qp_RNA_min = 2.23e-4;Qc_DNA = 9.41e-4;
Y_pc_DNA = 1/11.1;
A_pchl_Pho = 2.82e-2;
Qp_Other = 6.53e-4;
aP = Ap_RNA.*((Apho.*Achl)+Abio);
bP = Ap_RNA.*(Apho.*Bchl+Qc_Pro_Other).*Y_np_RNA+(A_pchl_Pho.*Achl);
cP =Qp_RNA_min+(Qc_DNA.*Y_pc_DNA)+(A_pchl_Pho.*Bchl)+Qp_Other;
mu_P=DSolver(aP,bP,cP);
Qp_Nsto = (aP.*mu_P.^2)+(bP.*mu_P)+cP;
Qp_sto = Qp-Qp_Nsto;
Qp_sto_max = 3.5e-2;
Qp_max = Qp_Nsto+Qp_sto_max;
Qp_exc = max(Qp-Qp_max,0);
Tp_exc = 4.17e-2;
kp = 3.65e0; Vpmax = 2.01e-3;
Vp = Vpmax.*(P./(P+kp));
Vp_eff = Vp-(Qp_exc./Tp_exc);
gQpm1 = Vp_eff-(mu_p.*Qp);
end