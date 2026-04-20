function [mu_N,mu_P,mu_Fe,mu,nc,pc,fec,np,nfe,...
         Vn,Vp,Vfe,Pchl] = alloc3(I,N,P,Fe,Pmax,kn,kp,kfe,vmaxn,vmaxp,vmaxfe,Qp_max)
% determining Qn
% constants related to Qn
Qc_prooth = 2.4e-1;
ync_pro = 1/4.49;
ynp_rna = 3.8/1;
ync_dna = 3.8/11.1;
ync_chl = 4/55;
Pchl=Pmax*(1-exp(-OT*I)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
Achl = (1+E)./Pchl;
Bchl = m./Pchl;
Qc_chl = Achl.*mu + Bchl;
Qc_propho = Apho.*Qc_chl; Qc_probio = Abio.*mu;
Qc_pro = Qc_propho + Qc_probio + Qc_prooth;
Qn = Qc_pro.*ync_pro + Qp_rna.*ynp_rna...
+ Qc_dna.*ync_dna + Qc_chl.*ync_chl;



end