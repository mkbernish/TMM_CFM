function [mu_N,mu_P,mu_Fe,mu,nc,pc,fec,np,nfe,...
         Vn,Vp,Vfe,Pchl,lim] = alloc(I,N,P,Fe,Pmax,kn,kp,kfe,vmaxn,vmaxp,vmaxfe,Qp_max,A_pho_Fe)

%%%%%%%%%%%%%%% parameter sets %%%%%%%%%%%%%%%%%%%%%%%%
%Pmax=0.00320513285659728;
OT=0.00863364097132997;
%m=3.79146798299876E-19   ;      %(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
m = 2.16e-2;

Molar_mass_DNA_AT_average=308.47  ;      %(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_DNA_CG_average=308.97   ;     %(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_AT_average=317.47    ;    %(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=324.97     ;   %(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
%E coli
CG_Ecoli=0.506  ;        %(dimensionless) from (http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016))
AT_Ecoli=1-CG_Ecoli;     %(dimensionless)
Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli ;    %(g mol-1) Molar mass of DNA unit
Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli  ;   %(g mol-1) Molar mass of RNA unit
RNA_DNA_mass_ratio=17.844/6.5239;  %(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli  ;  %(mol mol-1)
%Stoichiometric parameters for DNA and RNA
CG=0.563      ;
AT = 1-CG;AU = 1-CG;
C_CG_RNA = 19/2;N_CG_RNA = 8/2;
C_AU_RNA = 19/2;N_AU_RNA = 7/2;
%DNA
C_CG_DNA = 19/2;N_CG_DNA = 8/2;
C_AT_DNA = 20/2;N_AT_DNA = 7/2;
YnucacidN_P = N_CG_RNA*CG + N_AU_RNA*AU; %same for RNA and DNA
YnucacidP_N = 1/YnucacidN_P;
YdnaC_N = (C_CG_DNA*CG + C_AT_DNA*AT)/(N_CG_DNA*CG + N_AT_DNA*AT);YrnaC_N = (C_CG_RNA*CG + C_AU_RNA*AU)/(N_CG_RNA*CG + N_AU_RNA*AU);

DNAmb=2.1269     ;              %(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs (http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016))
Avogadro=6.022*10^23  ;         %(molecules mol-1) Avogadro constant
Pdna_const=DNAmb*2*10^6/Avogadro ;               %(molP cell-1) Constant part of DNA in phosphorus
%Prna_const=Pdna_const*RNA_DNA_molar_ratio   ;    %(molP cell-1) Constant part of RNA in phosphorus
Ndna_const=Pdna_const/YnucacidP_N  ;    %(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio ;  %(molN cell-1) Constatn part of RNA in phosphorus
%Ndna=Ndna_const  ;  %(molN cell-1) DNA in nitrogen (here assuming constant)
%%%%% mol N   mol C     mol N
%%%%% ----- X -----  -> -----
%%%%% mol C   cell       cell

%       mmol P
%   B : ------
%         m^3
E=0.7742;
Qc=1.00E-12/12 ;     %(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55;

%Conversion parameters================
CNprotein=4.49 ;  %(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanoC_N=2  ;                           %(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40   ;                        %(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
Feunit = 1/Qc;
YphotoFe_N=0.001636364  ;%(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)
Cnrna_variable=4.23e-3;
Cnbiosynth=0.271;
%Ynphoto_chl = 16;
Ypthylakoid_chl=0.0281633095303638 ;  
Nconst_protein = 0.24;
Prna_const = 2.23e-4;
Pconst_other = 6.53e-4;
Cdna_const = 9.41e-4;
%Cnrna_variable=6212.59249917364;    %(s) Constant for Variable part of RNA (193-26)
Pchl=Pmax*(1-exp(-OT*I)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
Achl = (1+E)./Pchl;
Bchl = m./Pchl;
%%%%%%% calculating iron growth rate (no storage) %%%%%%%%%%
% aff_Fe = 3.042e-17; % for e. hux - change depending on species?
Vfe = vmaxfe.*(Fe./(Fe+kfe));
%A_pho_Fe = 4e-3;
%Vfe=aff_Fe.*Fe    ;%(mol fe cell-1 s-1) iron uptake per cell
aFe = A_pho_Fe.*Achl;
bFe = A_pho_Fe.*Bchl;
cFe = -Vfe;
mu_Fe= real((-bFe+sqrt(bFe.^2-(4.*aFe.*(cFe))))./(2.*aFe));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% calculating nitrogen growth rate (no storage) %%%%%%%%%%
Ap_RNA = 4.23e-3;%/Arr;
Apho = 1.6e1;
Abio = 2.71e-1;
Y_np_RNA = 3.8/1;
Y_nc_Pro = 1/4.49;
Qc_Pro_Other = 2.4e-1;
Qp_RNA_min = 2.23e-4;
Qc_DNA = 9.41e-4;
Y_nc_DNA = 3.8/11.1;
Ynphoto_chl = 4/55;
Y_nc_chl = 4/55; 
%aff_N = 2.7e-20;
Vn = vmaxn.*(N./(N+kn));
%aN = Y_np_RNA.*Ap_RNA.*(Apho.*Achl+Abio);
% bN = Y_nc_chl.*Achl+Y_nc_Pro.*(Apho.*Achl+Abio)+...
%     Y_np_RNA.*Ap_RNA.*(Apho.*Bchl+Qc_Pro_Other);
% cN = Ynphoto_chl.*Bchl+Y_nc_Pro.*(Qc_Pro_Other+Apho.*Bchl)+...
%     Y_np_RNA.*Qp_RNA_min+Y_nc_DNA.*Qc_DNA;
aN = Ap_RNA.*((Apho.*Achl)+Abio).*Y_np_RNA;
bN = ((Apho.*Achl)+Abio).*Y_nc_Pro+(Achl.*Ynphoto_chl)+...
    Ap_RNA.*(Apho*Bchl+Qc_Pro_Other).*Y_np_RNA;
cN = Bchl.*Ynphoto_chl+(Apho*Bchl+Qc_Pro_Other).*...
    Y_nc_Pro+Qp_RNA_min*Y_np_RNA+Qc_DNA*Y_nc_DNA;
dN = -Vn;
mu_N=real(-bN./(3.*aN)...
-((2.^(1/3).*(-bN.^2+3.*aN.*cN))...
./(3.*aN.*(-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN+(4.*(-bN.^2+3.*aN.*cN).^3+...
(-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN).^2).^(1/2)).^(1/3)))...
+((-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN+(4.*(-bN.^2+3.*aN.*cN).^3+(-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN).^2).^(1/2)).^(1/3)...
./(3*2.^(1/3).*aN)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% calculating phosphorus growth rate (no storage) %%%%%%%%%%
Ap_RNA = 4.23e-3;

Abio = 2.71e-1;
Qc_Pro_Other = 2.4e-1;
Qp_RNA_min = 2.23e-4;
Qc_DNA = 9.41e-4;
Y_pc_DNA = 1/11.1;
Y_pchl_Pho = 2.82e-2;
Qp_Other = 6.53e-4;
Vp = vmaxp.*(P./(P+kp));
%aP = Ap_RNA.*(Apho.*Achl+Abio);
%bP = Achl.*Y_pchl_Pho+Ap_RNA.*(Apho.*Bchl+Qc_Pro_Other);
%cP = Qp_Other+(Bchl.*Y_pchl_Pho)+Qp_RNA_min+(Y_pc_DNA.*Qc_DNA);

aP = Ap_RNA.*((Apho.*Achl)+Abio);
bP = Ap_RNA.*(Apho.*Bchl+Qc_Pro_Other)+(Y_pchl_Pho.*Achl);
cP =Qp_RNA_min+(Qc_DNA.*Y_pc_DNA)+(Y_pchl_Pho.*Bchl)+Qp_Other;

dP = -Vp;
mu_P=real(-bP./(3.*aP)...
-((2.^(1/3).*(-bP.^2+3.*aP.*cP))...
./(3.*aP.*(-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP+(4.*(-bP.^2+3.*aP.*cP).^3+...
(-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP).^2).^(1/2)).^(1/3)))...
+((-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP+(4.*(-bP.^2+3.*aP.*cP).^3+(-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP).^2).^(1/2)).^(1/3)...
./(3*2.^(1/3).*aP)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solving carbon?
YproteinC_N = 4.49 ;        %(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
Phi_maint = 3.791468e-19;   
Cessential = 1.51786753491048e-15  ;        %(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
Ynphoto_chl = 3.56099164557551 ;  
Cnbiosynth = 4.34728279914354e-10 ; 
Cnrna_variable = 6212.59249917364 ;       %(s) Constant for Variable part of RNA (193-26)
Nconst_protein = 4.45336898828389e-15 ;  
vIMax = 0.00320513 ;
Pchl2=Pmax*(1-exp(-vIMax*I)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)

aC = (((1+E)*Qc*Ynphoto_chl)./Pchl2 + Cnbiosynth)*Cnrna_variable*YrnaC_N ; 
bC = (1+E)*Qc./Pchl2 ...
    + (((1+E)*Qc*Ynphoto_chl)./Pchl2 + Cnbiosynth)*YproteinC_N ...
    + (((1 + E)*Qc*Ypthylakoid_chl)./Pchl2)*YpgC_P ...
    + (Nconst_protein+(Phi_maint*Ynphoto_chl)./Pchl2)*Cnrna_variable*YrnaC_N ;

cC = Phi_maint./Pchl2 ...
    + (Nconst_protein+(Phi_maint*Ynphoto_chl)./Pchl2)*YproteinC_N ...
    + ((Phi_maint*Ypthylakoid_chl)./Pchl2)*YpgC_P ...
    + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential ...
    -Qc;
mu_C= real((-bC+sqrt(bC.^2-(4.*aC.*(cC))))./(2.*aC)).*(60*60*24);

%%%%%%% effective growth rate %%%%%%%%%%
mu2 = min(mu_N,mu_P);
mu = min(mu2,mu_Fe);
Qc_essential = aC.*mu.^2 + bC.*mu + cC ;

%%

[Nlim,~] = find(mu_N==mu);
[Plim,~] = find(mu_P==mu);
[Felim,~] = find(mu_Fe==mu);
lim = nan(size(mu_Fe));
lim(mu_Fe<mu_N & mu_Fe<mu_P) = 0;
lim(mu_N<mu_Fe & mu_N<mu_P) = 10;
lim(mu_P<mu_Fe & mu_P<mu_N) = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Chl       =Bchl + Achl.*mu;
Nchl=Chl*YchlN_C ;       % %(molN chl cell-1) Chlorophyll N concentration
Pthylakoid=Chl*Ypthylakoid_chl         ;% %(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
Nphoto_const = Bchl*Ynphoto_chl        ;    % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nphoto_D = Achl.*Ynphoto_chl            ;        % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Fephoto_D = Achl.*A_pho_Fe;
Fephoto_const = Bchl.*A_pho_Fe;
Nprotein_const = Nphoto_const + Nconst_protein;  % (molN cell-1) All the proteins in N (193-26)
Nprotein_D = Nphoto_D + Cnbiosynth  ;
Nprotein  = Nprotein_const + Nprotein_D.*mu ; %masked array, ~7.4 e-10     %(mol N cell-1) all the protein
Fephoto = Fephoto_D+Fephoto_const;
Ndna_variable=Ndna_const.*mu/1.2*(18.3-7.6)/7.6     ;%   %(molN cell-1) variable part of nitrogen in DNA (193-26) Increasing ratio based on Bremmer 1996
Nrna_variable=Nprotein .* mu .* Cnrna_variable;  %      %(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Punit = 1/Qc;
Pdna_variable=Ndna_variable*YnucacidP_N  .*Punit  ;% %(molP cell-1) variable part of phosphorus in DNA (193-26)
Prna_variable=Nrna_variable*YnucacidP_N ;%   %(molP cell-1) variable part of phosphorus in RNA (193-26)Pdna_variable = 1;

%%% doing it like in the paper


Qc_Chl = Achl.*mu+Bchl;
Qc_propho = Apho.*Qc_Chl;
Qc_probio = Abio.*mu;
Qc_Pro = Qc_propho + Qc_probio + Qc_Pro_Other;
Npro2 = Qc_Pro .* Y_nc_Pro;
Qp_RNA = Ap_RNA.*mu.*Qc_Pro+Qp_RNA_min;
Qn_RNA = Qp_RNA.*Y_np_RNA;
Qn_DNA = Qc_DNA.*Y_nc_DNA;
Qn_chl = Qc_Chl.*Y_nc_chl;
Qn=Npro2+Qn_RNA+Qn_DNA+Qn_chl ;%  
Pdna_const_plot = Pdna_const.*Punit;
Qp1=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const_plot   ;% 
Qp_DNA = Qc_DNA.*Y_pc_DNA;
Apho_pchl = Y_pchl_Pho;
Qp_thy = Apho_pchl.*Qc_Chl;
Qp = Qp_RNA + Qp_DNA + Qp_thy + Qp_Other;
Qfe1 = Fephoto;
Qfe_pho = A_pho_Fe.*Qc_Chl;
Qfe = Qfe_pho;

Qntot = Vn./mu;
Qptot = Vp./mu;
Qfetot = Vfe./mu;

Qn_sto_max = 3.5e-2;
Qn_sto = max(0,Qntot-Qn);
idx = Qn_sto>Qn_sto_max; Qn_sto(idx) = Qn_sto_max;
Qn3 = Qn+Qn_sto;

Qp_sto = max(0,Qptot-Qp);
idx = Qptot>Qp_max; 
Qptot(idx) = Qp_max;
%idx = Qp> Qp_max;
%Qp(idx) = Qp_max;
Qfe_max = 2.44e-4;
Qp3 = Qptot;
idx = Qfetot>Qfe_max;
Qfetot(idx) = Qfe_max;
Qfe3 = Qfetot;
%%%%%%% non-storage molecules %%%%%%%%%%
Qn_Nsto =  aN.*mu.^2+bN.*mu+cN;
Qp_Nsto =  aP.*mu.^2+bP.*mu+cP;
 Qfe_Nsto =  aFe.*mu+bFe; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qn(Nlim) =  aN(Nlim).*mu(Nlim).^2+bN(Nlim).*mu(Nlim)+cN(Nlim);
%Qp(Plim) =  aP(Plim).*mu(Plim).^2+bP(Plim).*mu(Plim)+cP(Plim);
%Qfe(Felim) = aFe(Felim).*mu(Felim)+(bFe(Felim)); 

%idx = Qn> (Qn_sto_max+Qn_Nsto);
%Qn(idx) = Qn_sto_max+Qn_Nsto(idx);


% idx = Qfe>Qfe_max;
% Qfe(idx) = Qfe_max;


%%%%%%% Total quota? %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% storage molecules %%%%%%%%%%

Qn_sto =  Qn-Qn_Nsto;
idx = Qn_Nsto> Qn; Qn_sto(idx) = 0;
idx = Qn_sto<1e-15; Qn_sto(idx) = 0;

Qn_max = Qn_Nsto+Qn_sto_max;
stoNdiff = zeros(size(Qn));
if any(Qn_sto>Qn_sto_max)
    idx = Qn_sto>Qn_sto_max; 
    stoNdiff(idx) = Qn_sto(idx)-Qn_sto_max;
    Qn_sto(idx)=Qn_sto_max;
end
Qn_exc = max(Qn-Qn_max+stoNdiff,0);

idx = Qn_Nsto> Qn; Qn_sto(idx) = 0;
idx = Qn_sto<1e-15; Qn_sto(idx) = 0;
Qn_sto_max = 3.5e-2;
Qp_sto =  Qp-Qp_Nsto;
idx = Qp_Nsto>Qp;Qp_sto(idx) = 0; % when the p allocated to growth is greater than
% the total P allocation (for some reason), sets P storage to 0
idx = Qp_sto<1e-16; Qp_sto(idx) = 0;
idx = Qp>Qp_max;Qp(idx)=Qp_max;
Qfe_sto =  Qfe-Qfe_Nsto;
idx = Qfe_Nsto> Qfe;Qfe_sto(idx) = 0;
idx = Qfe_sto<1e-16; Qfe_sto(idx) = 0;
%idx = Qfe>Qfe_max;
%Qfe(idx)=Qfe_max;

Qn_exc = max(Qn-Qn_max+stoNdiff,0);
Tn_exc = 4.17e-2;
Vn = vmaxn.*(N./(N+kn));
Vn_eff = Vn-(Qn_exc./Tn_exc);
idx = N<=0; Vn_eff(idx) = 0;
%gQn = Vn_eff-((mu).*Qn);
% Qn2 = Vn./mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np=Qn3./Qp3 ;%   
%np=Qn./Qp; 
nc=Qn3; 
pc=Qp3; fec=Qfe3;nfe=Qn3./Qfe3;
%nc = Qn_Nsto;
end