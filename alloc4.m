function [out] = alloc4(I,N,P,Fe,Pmax,...
    kn,kp,kfe,vmaxn,vmaxp,vmaxfe,Qp_max,A_pho_Fe)

%%%%%%%%%%%%%%% parameter sets %%%%%%%%%%%%%%%%%%%%%%%%
OT=0.00863364097132997;

m=3.79146798299876E-19   ;      %(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
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
Ndna_const=Pdna_const/YnucacidP_N  ;    %(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio ;  %(molN cell-1) Constatn part of RNA in phosphorus
Ndna=Ndna_const  ;  %(molN cell-1) DNA in nitrogen (here assuming constant)

E=0.7742;
Qc=1.00*10^(-12)/12 ;     %(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
vIMax = Pmax;
%Conversion parameters================
Pchl=vIMax.*(1-exp(-OT*I)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
Achl = (1+E)./Pchl; %Chl_D
Bchl = m./Pchl; %Chl_const
%%

%%%%% put things from gabbys slide here
%Cnrna_variable=6212.59249917364;    %(s) Constant for Variable part of RNA (193-26)
%Cnbiosynth=4.34728279914354E-10 ;      %(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
%Ynphoto_chl = 3.56099164557551; %((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
%Nconst_protein = 4.45336898828389E-15 ;  %(molN cell-1) Constant protein pool in nitrogen (193-25)
%Prna_const=Pdna_const*RNA_DNA_molar_ratio   ; %Qp RNA Min    %(molP cell-1) Constant part of RNA in phosphorus
%Pconst_other=5.44534485638617E-17  ;             %(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
%Cdna_const = Ndna_const.*YdnaC_N;

%% kei's paper
Ynphoto_chl = 3.56099164557551 ;         %((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
Cnbiosynth = 4.34728279914354e-10 ;       %(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Nconst_protein = 4.45336898828389e-15 ;     %(molN cell-1) Constant protein pool in nitrogen (193-25)
Cnrna_variable = 6212.59249917364 ;       %(s) Constant for Variable part of RNA (193-26)
YchlN_C = 4/55 ;
% P vars
Ypthylakoid_chl = 0.0281633095303638 ;       %((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other = 5.44534485638617e-17 ;              %(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
YproteinC_N = 4.49 ;        %(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanophycinC_N = 2  ;      % (molC molN) C/N molar ratio of cyanophycin (N storage molecule) - Cyanophycin is a carbon/nitrogen storage polymer found in cyanobacteria and in a few heterotrophic bacteria.
YpgC_P = 40  ;              % (molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (phospholipid in thylakoid membrane) (assumes C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
% C vars
Cessential = 1.51786753491048e-15  ;        %(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
% Fe Vars
YphotoFe_N = 0.001636364  ;   %(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)
vIMax = Pmax;
ksatPAR    = 0.008633641;   % (uEin m^-2 s^-1)^-1 light dependence of chlorophyll carbon fixation
Pchl = vIMax.*(1 - exp(-ksatPAR*I));
%%%%%%% calculating phosphorus growth rate (no storage) %%%%%%%%%%
%Vp2 = vmaxp.*(P./(P+kp)).*Qc./86400;

%% megan P - calculating mu
%aPO4 = 0.01*Qc/86400  ;%(mol N cell-1 s-1 / umolN L-1) affinity of NO3 (initial value from d007_06_00)
%affP = vmaxp/kp;
%affinityPO4 = affP.*Qc./86400;
%Vp2 = aPO4.*P ;
Vp2 = vmaxp.*(P./(P+kp)).*Qc./86400;
aP = (((1+E)*Qc*Ynphoto_chl)./Pchl + Cnbiosynth )*Cnrna_variable*YnucacidP_N ;
bP = ((1+E)*Qc*Ypthylakoid_chl)./Pchl ...
    + (Nconst_protein +(m*Ynphoto_chl)./Pchl)*Cnrna_variable*YnucacidP_N ;
cP = (m*Ypthylakoid_chl)./Pchl + Nrna_const*YnucacidP_N + Ndna*YnucacidP_N +Pconst_other ;
dP = -Vp2 ;
mu_P=real(-bP./(3.*aP)...
    -((2.^(1/3).*(-bP.^2+3.*aP.*cP))...
    ./(3.*aP.*(-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP+(4.*(-bP.^2+3.*aP.*cP).^3+...
    (-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP).^2).^(1/2)).^(1/3)))...
    +((-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP+(4.*(-bP.^2+3.*aP.*cP).^3+(-2.*bP.^3+9.*aP.*bP.*cP-27.*aP.^2.*dP).^2).^(1/2)).^(1/3)...
    ./(3*2.^(1/3).*aP)));

%% megan N
%affN = vmaxn/kn;
%affinityNO3 = affN.*Qc./86400;
%Vn = affinityNO3 .* N ;
Vn = vmaxn.*(N./(N+kn)).*Qc./86400;
aN = ( ((1+E)*Qc*Ynphoto_chl)./Pchl + Cnbiosynth )*Cnrna_variable ;
bN= ( ((1+E)*Qc*Ynphoto_chl)./Pchl +Cnbiosynth ) ...
    + (Nconst_protein + ( m*Ynphoto_chl)./Pchl )*Cnrna_variable ...
    + ((1+E)*Qc*YchlN_C)./Pchl ;
cN= (Nconst_protein + (m*Ynphoto_chl)./Pchl) +Nrna_const + Ndna + (m*YchlN_C)./Pchl ;
dN= -Vn ;
mu_N=real(-bN./(3.*aN)...
    -((2.^(1/3).*(-bN.^2+3.*aN.*cN))...
    ./(3.*aN.*(-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN+(4.*(-bN.^2+3.*aN.*cN).^3+...
    (-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN).^2).^(1/2)).^(1/3)))...
    +((-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN+(4.*(-bN.^2+3.*aN.*cN).^3+(-2.*bN.^3+9.*aN.*bN.*cN-27.*aN.^2.*dN).^2).^(1/2)).^(1/3)...
    ./(3*2.^(1/3).*aN)));

%% megan Fe
%affFe = vmaxfe/kfe;
%affinityFe = affFe.*Qc./86400;
%VFe = affinityFe .* Fe ;
VFe = vmaxfe.*(Fe./(Fe+kfe)).*Qc./86400;
YFe_photo1 = Qc*Ynphoto_chl*YphotoFe_N;
YFe_photo = Qc*A_pho_Fe;

%aFe = A_pho_Fe.*Achl;
aFe = YFe_photo* ((1+E)./Pchl ) ;
bFe = YFe_photo* (m./Pchl) ;
cFe= -VFe ;
mu_Fe= real((-bFe+sqrt(bFe.^2-(4.*aFe.*(cFe))))./(2.*aFe));

%% megan C
% aC = (((1+E)*Qc*Ynphoto_chl)./Pchl + Cnbiosynth)*Cnrna_variable*YrnaC_N ;
% bC = (1+E)*Qc./Pchl ...
%     + (((1+E)*Qc*Ynphoto_chl)./Pchl + Cnbiosynth)*YproteinC_N ...
%     + (((1 + E)*Qc*Ypthylakoid_chl)./Pchl)*YpgC_P ...
%     + (Nconst_protein+(m*Ynphoto_chl)./Pchl)*Cnrna_variable*YrnaC_N ;
% 
% cC = m./Pchl ...
%     + (Nconst_protein+(m*Ynphoto_chl)./Pchl)*YproteinC_N ...
%     + ((m*Ypthylakoid_chl)./Pchl)*YpgC_P ...
%     + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential ...
%     -Qc;
% 
% mu_C= real((-bC+sqrt(bC.^2-(4.*aC.*(cC))))./(2.*aC));

%% effective growth rate
mu = min([mu_N, mu_P, mu_Fe],[],2) ;
%mu = min([mu_N; mu_P; mu_Fe; mu_C],[],1) ;
limType = nan(size(N));
for ii = 1:length(N)
    if mu(ii) == mu_N(ii)
        limType(ii) = 1 ;    % N limitation
    elseif mu(ii) == mu_P(ii)
        limType(ii) = 2 ;    % P limitation
    elseif mu(ii) == mu_Fe(ii)
        limType(ii) = 3 ;    % Iron Limitation
    end
end

%% allocation - non-storage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chlorphyll related
Chl=((1+E)*mu*Qc+m)./Pchl   ;    %(molC chl cell-1) cN[i]hlrophyll concentration (193-25)
Chl_const = m./Pchl       ;                       % (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
Chl_D = (1 + E)*Qc./Pchl      ;                     % (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
Cchl_const = m./Pchl       ;                       % (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
Cchl_D = (1 + E)*Qc./Pchl      ;                     % (molC chl cell-1) cN[i]hlrophyll concentration (193-25)

% Nitrogen related
Nchl_D = Chl_D*YchlN_C   ;                       % (molN chl cell-1) Chlorophyll N concentration
Nchl_const = Chl_const*YchlN_C    ;              % (molN chl cell-1) Chlorophyll N concentration
Nprotphoto_D = Chl_D*Ynphoto_chl     ;               % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nprotphoto_const = Chl_const*Ynphoto_chl   ;         % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nprotbiosynth_D = Cnbiosynth             ;           % (molN cell-1) various part of biosynthesis related protein in N (193-37)
Nprotein_D = Nprotphoto_D + Nprotbiosynth_D       ;      % (molN cell-1) All the proteins in N (193-26)
Nprotein_const = Nprotphoto_const + Nconst_protein ; % (molN cell-1) All the proteins in N (193-26)
Nrna_D = Nprotein_const*Cnrna_variable   ;       % (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Nrna_D2 = Nprotein_D*Cnrna_variable     ;        % (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Nchl           = Nchl_const     + Nchl_D.*mu         ;      %(molN chl cell-1) Chlorophyll N concentration
Nprot_photo    = Nprotphoto_const   + Nprotphoto_D.*mu       ; %(molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nprot_biosynth =                  Nprotbiosynth_D.*mu  ;
Nprotein       = Nprotein_const + Nprotein_D.*mu ;
Nrna           = Nrna_const     + Nrna_D.*mu + Nrna_D2.*mu.*mu ;
Pdna=Ndna*YnucacidP_N ;  %(mol P cell-1) DNA in phosphorus
Pthylakoid=Chl*Ypthylakoid_chl  ;        %(molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
Prna=Nrna*YnucacidP_N ;    %(molP cell-1) Phosphorus in RNA
Fephoto_D =Chl_D.*A_pho_Fe.*mu;
Fephoto_const = Chl_const.*A_pho_Fe;
Fephoto = Fephoto_D+Fephoto_const;
%Pconst_other = Pconst_other ;
%=============================
%preparing for Qc computation
%=============================
%Constant carbon parameters------------------
Cconst_protein = Nconst_protein*YproteinC_N ; %(molC cell-1) carbon in other protein assumed constant (195-16)
Cdna_const = Ndna_const*YdnaC_N  ;    %(molC cell-1) carbon in constant part of DNA (195-16)

Crna = Nrna*YrnaC_N ;      %(molC cell-1) carbon in RNA (195-16)
Cchl = Cchl_const + Cchl_D.*mu   ;                %(molC cell-1) carbon in chlorophyll (195-16)
Cprot_photo = Nprot_photo*YproteinC_N   ;   %(molC cell-1) carbon in photosystem protein (195-16)
Cprot_biosynth = Nprot_biosynth*YproteinC_N ;   %(molC cell-1) carbon in biosynthesis protein (195-16)
CthylakoidPG = Pthylakoid*YpgC_P  ;       %(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
Qc_essential = Crna + Cdna_const + Cchl + Cprot_photo + Cprot_biosynth +Cconst_protein + CthylakoidPG + Cessential ;

 A=((1+E)*Qc*Ynphoto_chl)./Pchl+Cnbiosynth;
    B=Nconst_protein+(m*Ynphoto_chl)./Pchl;
    G=((1+E)*Qc*YchlN_C)./Pchl;
    H=(m*YchlN_C)./Pchl;
    I=A*Cnrna_variable;
    J=A+B*Cnrna_variable+G;
    K=B+Nrna_const+Ndna+H;
    L = ((1 + E)*Qc*Ypthylakoid_chl)./Pchl;
    M = (m*Ypthylakoid_chl)./Pchl;
    N = A*Cnrna_variable*YnucacidP_N;
    O = L + B*Cnrna_variable*YnucacidP_N;
    P = M + Nrna_const*YnucacidP_N + Ndna*YnucacidP_N + Pconst_other;

Qc_D2 = A*Cnrna_variable*YrnaC_N;
Qc_D = (1+E)*Qc./Pchl + A*YproteinC_N + L*YpgC_P + B*Cnrna_variable*YrnaC_N;
Qc_const = m./Pchl + B*YproteinC_N + M*YpgC_P + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential;

%mu2 = mu;
idx = Qc_essential >= Qc  ;%conditions where Mu max applies
mu2(idx) = 2*(Qc - Qc_const(idx))./(Qc_D(idx) + sqrt(Qc_D(idx).*Qc_D(idx) + 4*Qc_D2(idx).*(Qc - Qc_const(idx))));
Cother = max(0,Qc - Qc_essential);

Cnstore = zeros(size(N));
Nstor_max = 2.91679384515998E-15  ;  % this value is pretty close to half of the minimum functional N quota
% minimum N quota = total non-storage N quota
%Qn_nostor 
%Qn_nostor = Nprot_photo + Nprot_biosynth + Nconst_prot + Nchl + Ndna*makevec + Nrna ;
Qn_nsto = Nprot_photo + Nprot_biosynth + Nprotein + Nchl + Ndna*ones(size(Nprot_photo)) + Nrna ;
Qn_max = Qn_nsto+Nstor_max;

%Qn_min = Nprot_photo + Nprot_biosynth + Nprotein + Nchl + Ndna*ones(size(Nprot_photo)) + Nrna ;
Nstor = zeros(size(N)); % initialize vector of zeros
indx = find(limType~=1); % find not N limited locations

Nstor(indx) = max(0,Vn(indx)./mu(indx) - Qn_nsto(indx)); % calculate max N uptake for min mu
Nstor(Nstor > Nstor_max) = Nstor_max; % reset to maximum value
% this reset is non-linear. replace with a logistic function or similar to
% optimize.

YnstorC_N = YcyanophycinC_N ; % (cyanophycin)
Cnstor = Nstor.*YnstorC_N ;    % (molC cell-1) carbon in nitrogen storage (cyanophycin) 

% --- recalculate Cother ---
Cother = max(0,Qc - Qc_essential - Cnstor);
Qp_nsto=Pconst_other+Pthylakoid+Prna+Pdna  ;    %(molP cell-1) total phosphorus in the cell without storage
%Pstor_max2 = Nstor_max * (Vp2./Vn);
%Qp_max2 = Vp2./mu;
Qp_max1 = Qp_max.*Qc.*ones(size(P)); 
%Qp_max1 = Qp_nsto.*4;
%idx1 = find(limType==1);
%Qp_max1(idx1) = Qp_nsto(idx1).*2;
%
Pstor_max = Qp_max1-Qp_nsto;
%Pstor_max = Nstor_max * (Vp2./Vn);
Pstor = zeros(size(P)); % initialize vector of zeros
indx = find(limType~=2); % find not P limited locations
Pstor(indx) = Vp2(indx)./mu(indx) - Qp_nsto(indx); % calculate max P uptake for min mu; 
% max Qp = VP/mu = max potential P accumulation given the nutrient uptake rate 
% and growth rate
indx_maxPstor = find(Pstor > Pstor_max); 
Pstor(indx_maxPstor) = Pstor_max(indx_maxPstor); % reset to maximum value
%Qfe_Nsto = (A_pho_Fe.*Achl.*mu)+(A_pho_Fe.*Bchl);
Qfe_Nsto = Fephoto;
idx = find(limType~=3); % find not P limited locations
Qfe_sto = zeros(size(P)); 
Qfe_sto(idx) = (VFe(idx)./mu(idx))-Qfe_Nsto(idx);
Qfe_max = 2.44e-4;
Qfe_sto_max = Qfe_max - Qfe_Nsto;
idx = Qfe_sto>Qfe_sto_max;
Qfe_sto(idx)=Qfe_sto_max(idx);
%Qfe_nsto = 
out.mu = mu ;
out.limType = limType ; 
out.Qc = Qc ;
out.Cprot_photo = Cprot_photo ;
out.Cprot_biosynth = Cprot_biosynth ;
out.Cconst_prot = Cconst_protein ;
out.Cchl = Cchl ;
out.Crna = Crna ;
out.Cdna_const = Cdna_const ;
out.Cother = Cother ;
out.Cessential = Cessential ;
out.CthylakoidPG = CthylakoidPG ;
out.Cnstore = Cnstore ;
out.Nchl = Nchl ;
out.Nprot_photo = Nprot_photo ;
out.Nprot_biosynth = Nprot_biosynth ;
out.Nconst_prot = Nconst_protein ;
out.Nprotein = Nprotein ;
out.Nrna = Nrna ;
out.Ndna = Ndna ;   
out.Nstor = Nstor ;

out.Pthylakoid = Pthylakoid ;
out.Pdna = Pdna ;
out.Prna = Prna ;
out.Pconst_other = Pconst_other;
out.Pstor = Pstor ;
out.Qfe_Nsto = Qfe_Nsto;
out.Qfe_sto = Qfe_sto;
makevec = ones(size(N));
out.Qc_essential = Qc_essential;
%out.Qc_functional = Qc_essential; 
out.Qn_nostor = out.Nprot_photo + out.Nprot_biosynth + out.Nconst_prot + out.Nchl + out.Ndna*makevec + out.Nrna ;
out.Qp_nostor = out.Pthylakoid + out.Pdna*makevec + out.Prna + out.Pconst_other*makevec;
out.Qn_total = out.Qn_nostor + out.Nstor ;
out.Qp_total = out.Qp_nostor + out.Pstor ;
out.Qfe_total = out.Qfe_Nsto+out.Qfe_sto;
%idx = out.Qp_total>Qp_max1;
%out.Qp_total(idx) = Qp_max1(idx);
out.CP = out.Qc./out.Qp_total;
out.NP = out.Qn_total./out.Qp_total ;
out.CN = out.Qc./out.Qn_total ; 
out.PC = out.Qp_total./out.Qc;
out.NC = out.Qn_total./out.Qc; 
out.FeC = out.Qfe_total./out.Qc;

out.fracPhoto = (out.Cprot_photo + out.Cchl + out.CthylakoidPG)./Qc;
out.fracBiosynth = (out.Cprot_biosynth + out.Crna)./Qc;
out.fracStruct_const = (out.Cdna_const.*makevec + out.Cessential.*makevec + out.Cconst_prot.*makevec)./Qc ; % fixed fraction of cell C in essential structure pools
out.fracStor = (out.Cother + out.Cnstore)./Qc ; 
out.Vp = Vp2;
out.Vn = Vn;
out.mu_N = mu_N;
out.mu_P = mu_P;
%out.mu_C = mu_C;
out.mu_Fe = mu_Fe;
out.Pchl = Pchl;
out.Pstor_max = Pstor_max;
%%%%%%% non-storage molecules %%%%%%%%%%

end