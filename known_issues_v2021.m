%% Known issues v2.2021
clear all 
close all
clc

% Load v2.2020 merged master file
path='C:\Users\nlange\Downloads';
cd(path)
load('GLODAPv2.2020_Merged_Master_File.mat')

% Get parameter names
params=who;
params([3,4,9,26,27,28,46,47,48,49,50,51,92,105,106,107,108])=[];
params(8)=[];

%% Rows with missing Temperature to NaN
% 6532 samples 
% 130 different cruises 

ind=find(isnan(G2temperature) & (~isnan(G2salinity) | ~isnan(G2oxygen) | ~isnan(G2phosphate) | ~isnan(G2nitrate) | ~isnan(G2silicate) | ~isnan(G2cfc11) | ~isnan(G2cfc12) | ~isnan(G2cfc113)));
for j=1:length(params)
    if strcmp(params{j}(end),'f')==1   
        str=strcat(params{j},'(ind)=9;');
        eval(str)
        clear str
    elseif strcmp(params{j}(end-1:end),'qc')==1   
        str=strcat(params{j},'(ind)=',params{j},'(ind);');
        eval(str)
        clear str
    else   
        str=strcat(params{j},'(ind)=NaN;');
        eval(str)
        clear str
    end
end
clear ind


%% Fix "easy" bugs, usually meaning to delete bad samples
% Salinity
cruise=[21,245,268,362,459,1039];
station=[15,105,287,206,27,126];
bottle=[37,48,39,2,4,12];
cast=[1,2,1,1,1,1];

for i=1:length(cruise)
    ind=find(G2cruise==cruise(i) & G2station==station(i) & G2cast==cast(i) & G2bottle==bottle(i));
    if isempty(ind) | length(ind)>1
        disp(strcat('Warning for cruise:',num2str(cruise(i))))
    else
        G2salinity(ind)=NaN;
        G2salinityf(ind)=9;
    end
    clear ind
end
clear cruise i station bottle cast 


% Nutrients
cruise=[21,255];
station=[107,1005];
bottle=[15,23];


for i=1:length(cruise)
    ind=find(G2cruise==cruise(i) & G2station==station(i) & G2bottle==bottle(i));
    if isempty(ind) | length(ind)>1
        disp(strcat('Warning for cruise:',num2str(cruise(i))))
    else
        G2phosphate(ind)=NaN;
        G2phosphatef(ind)=9;
        G2silicate(ind)=NaN;
        G2silicatef(ind)=9;
        G2nitrate(ind)=NaN;
        G2nitratef(ind)=9;
    end
    clear ind
end

clear cruise i station bottle cast 


% Oxygen
cruise=[235];
station=[18];
bottle=[3];
cast=[2];

for i=1:length(cruise)
    ind=find(G2cruise==cruise(i) & G2station==station(i) & G2cast==cast(i) & G2bottle==bottle(i));
    if isempty(ind) | length(ind)>1
        disp(strcat('Warning for cruise:',num2str(cruise(i))))
    else
        G2oxygen(ind)=NaN;
        G2oxygenf(ind)=9;
    end
    clear ind
end

clear cruise i station bottle cast 

cruise=[279,279,279,279,279,279,279];
station=[1:7];

for i=1:length(cruise)
    ind=find(G2cruise==cruise(i) & G2station==station(i));
    if isempty(ind) | length(ind)>24
        disp(strcat('Warning for station:',num2str(station(i))))
    else
        G2oxygen(ind)=NaN;
        G2oxygenf(ind)=9;
    end
    clear ind
end

clear cruise i station bottle cast 

% TA
cruise=[238,345,345,345];
station=[32,94,98,118];
bottle=[1,12023,12020,12025];

for i=1:length(cruise)
    ind=find(G2cruise==cruise(i) & G2station==station(i) & G2bottle==bottle(i));
    if isempty(ind) | length(ind)>1
        disp(strcat('Warning for cruise:',num2str(station(i))))
    else
        G2talk(ind)=NaN;
        G2talkf(ind)=9;
    end
    clear ind
end

clear cruise i station bottle cast 

cruise=656;

ind=find(G2cruise==cruise);
G2talk(ind)=NaN;
G2talkf(ind)=9;
G2phts25p0(ind)=NaN;
G2phts25p0f(ind)=9;
G2phtsinsitu(ind)=NaN;
G2phtsinsituf(ind)=9;
G2fco2(ind)=NaN;
G2fco2f(ind)=9;

clear cruise i station bottle cast 


% Temperature
cruise=[268];
station=[226];
bottle=[94];
cast=[8];

for i=1:length(cruise)
    ind=find(G2cruise==cruise(i) & G2station==station(i) & G2cast==cast(i) & G2bottle==bottle(i));
    if isempty(ind) | length(ind)>1
        disp(strcat('Warning for cruise:',num2str(cruise(i))))
    else
        for j=1:length(params)
            if strcmp(params{j}(end),'f')==1
                str=strcat(params{j},'(ind)=9;');
            elseif strcmp(params{j}(end-1:end),'qc')==1
                str=strcat(params{j},'(ind)=',params{j},'(ind);');
            else
                str=strcat(params{j},'(ind)=NaN;');
            end
            eval(str)
            clear str
        end                
    end
    clear ind
end

clear cruise i station bottle cast params


%% Solve "bad" qc flags
% This includes qc-flags=1 and only NaN-values (only present up to cruiseno
% 724, i.e v2-only) and qc-flags=1 and only calculated values
cruises=unique(G2cruise,'stable');
params={'G2salinity','G2tco2','G2talk','G2phts25p0','G2nitrate','G2phosphate','G2silicate','G2oxygen','G2cfc12','G2cfc11','G2cfc113','G2ccl4','G2c13'};
for i=1:length(cruises)
    ind=find(G2cruise==cruises(i));
    for j=1:length(params)
        if j==4
            str=strcat('varf=',params{j},'f(ind);');
            eval(str)
            clear str
            str=strcat('varqc=',params{j}(1:6),'qc(ind);');
            eval(str)
            clear str
        else
            str=strcat('varf=',params{j},'f(ind);');
            eval(str)
            clear str
            str=strcat('varqc=',params{j},'qc(ind);');
            eval(str)
            clear str
        end
        % If qc-flag=1 but only calculated or NaN values present for the
        % cruise change qc-flag to 0!
        if max(varqc)==1 & max(ismember(unique(varf),2))==0
            if j==4
                str=strcat(params{j}(1:6),'qc(ind)=0;');
                eval(str)
                clear str
            else
                str=strcat(params{j},'qc(ind)=0;');
                eval(str)
                clear str
            end
        end
        clear var* 
    end
    clear ind
end
clear i j params

%% Give -777 adjustment a qc-flag of 1 "back"
% If adjustment table shows a -777 for a parameter fo a cruise, a qc-flag=1 should be assigned... 
% Load adjustment table first
T=readtable('GLODAPv2_adjustment_table.csv');
for i=1:size(T,1)
    holder=table2array(T(i,8:29));
    bflag=find(holder==-777);
    if isempty(bflag)
        continue
    end
    no=find(strcmp(expocode,T.cruise_expocode(i))==1);
    cruiseno=expocodeno(no);
    ind=find(G2cruise==cruiseno);
    
    if max(ismember(bflag,[1]))==1
        G2salinityqc(ind)=1;
    end
    if max(ismember(bflag,[6]))==1
        G2tco2qc(ind)=1;
    end
    if max(ismember(bflag,[7]))==1
        G2talkqc(ind)=1;
    end
    if max(ismember(bflag,[8]))==1
        G2phtsqc(ind)=1;
    end
    if max(ismember(bflag,[10]))==1
        G2nitrateqc(ind)=1;
    end
    if max(ismember(bflag,[11]))==1
        G2phosphateqc(ind)=1;
    end
    if max(ismember(bflag,[12]))==1
        G2silicateqc(ind)=1;
    end
    if max(ismember(bflag,[13]))==1
        G2oxygenqc(ind)=1;
    end
    if max(ismember(bflag,[18]))==1
        G2cfc12qc(ind)=1;
    end
    if max(ismember(bflag,[19]))==1
        G2cfc11qc(ind)=1;
    end
    if max(ismember(bflag,[20]))==1
        G2cfc113qc(ind)=1;
    end
    if max(ismember(bflag,[21]))==1
        G2ccl4qc(ind)=1;
    end
    if max(ismember(bflag,[22]))==1
        G2c13qc(ind)=1;
    end
    clear ind no cruiseno holder bflag
end
clear bflag holder i T

    
%% Assign 20°C temperature to calculated fco2 data (where missing)
ind=find(~isnan(G2fco2) & isnan(G2fco2temp)); % When checked, all these values have a flag=0
G2fco2temp(ind)=20;
clear ind


%% 2nd QC updates for existing cruises
% Give 316N19950124 a 5umol/kg downward adjustment for TA
% ind=find(G2cruise==250);
% G2talk(ind)=G2talk(ind)-5;
% 
% %Recalculate pH-tot with new TA values
% [DATA,~,~]=CO2SYS(G2talk(ind),G2tco2(ind),1,2,G2salinity(ind),G2temperature(ind),25,G2pressure(ind),0,G2silicate(ind),G2phosphate(ind),1,10,1);
% G2phts25p0(ind)=DATA(:,37);
% clear DATA
% 
% %Recalculate pH-inSitu with new TA values
% [DATA,~,~]=CO2SYS(G2talk(ind),G2tco2(ind),1,2,G2salinity(ind),G2temperature(ind),G2temperature(ind),G2pressure(ind),G2pressure(ind),G2silicate(ind),G2phosphate(ind),1,10,1);
% G2phtsinsitu(ind)=DATA(:,37);
% 
% clear ind cruises path DATA
% 
% 
% %Give 316N19950830 a 5umol/kg downward adjustment for TA
% ind=find(G2cruise==255 & G2station>=962);
% G2talk(ind)=G2talk(ind)-5;
% 
% %Recalculate pH-tot with new TA values
% [DATA,~,~]=CO2SYS(G2talk(ind),G2tco2(ind),1,2,G2salinity(ind),G2temperature(ind),25,G2pressure(ind),0,G2silicate(ind),G2phosphate(ind),1,10,1);
% G2phtsinsitu(ind)=DATA(:,37);
% clear DATA
% 
% %Recalculate pH-inSitu with new TA values
% [DATA,~,~]=CO2SYS(G2talk(ind),G2tco2(ind),1,2,G2salinity(ind),G2temperature(ind),G2temperature(ind),G2pressure(ind),G2pressure(ind),G2silicate(ind),G2phosphate(ind),1,10,1);
% G2phts25p0(ind)=DATA(:,37);
% 
% clear ind cruises path DATA


% Give 18SN19950803 a 8% downward adjustment for PO4
ind=find(G2cruise==198);
G2phosphate(ind)=G2phosphate(ind).*0.92;
clear ind

% Give 49NZ20020822 a 6% upward adjustment for PO4
ind=find(G2cruise==480);
G2phosphate(ind)=G2phosphate(ind).*1.06;
clear ind

% % Give 49NZ20130106 a 2% upward adjustment for SiOH4
% ind=find(G2cruise==1051);
% G2silicate(ind)=G2silicate(ind).*1.02;
% clear ind

    
%% Update 33RO19980123 TA for Station 106
cd 'D:\GLODAPv2.2021\Updates\33RO19980123'
A=load('33RO19980123.mat');

ind=find(G2cruise==341);
A.ALKALI(A.ALKALI_FLAG_W~=2 & A.ALKALI_FLAG_W~=6)=NaN;
A.ALKALI_FLAG_W(A.ALKALI_FLAG_W~=2 & A.ALKALI_FLAG_W~=6)=9;

inds=find(A.STNNBR==106);
indsg=find(G2station(ind)==106);

G2talk(ind(indsg))=A.ALKALI(inds);
G2talkf(ind(indsg))=A.ALKALI_FLAG_W(inds);

clear ind inds indsg A

%% Update carbon data and TDN for 33RR20160208 (received new TDN, pH and TA data)
cd 'D:\GLODAPv2.2021\Updates\33RR20160208'
A=load('33RR20160208.mat');

% Sort the same way GLODAP data is sorted, i.e. first using station and
% then pressure
A=struct2table(A);
A=sortrows(A, [68 16]);

% Change flags for TDN
A.TDN(A.TDN_FLAG_W~=2 & A.TDN_FLAG_W~=6)=NaN;
A.TDN_FLAG_W(A.TDN_FLAG_W~=2 & A.TDN_FLAG_W~=6)=9;

% Convert pH for stations 25 and 26 from 9°C to 25°C
indph=find(A.STNNBR==25 | A.STNNBR==26);

% Calculate missing PO4 and SIOH4 values
SILCAT_ca=zeros(length(A.SILCAT),1);
PHSPHT_ca=zeros(length(A.PHSPHT),1);
A.YEAR=A.SALNTY;
A.YEAR(:)=2016;
data=[A.YEAR,A.LATITUDE,A.LONGITUDE,A.CTDPRS,A.CTDTMP,A.SALNTY,A.OXYGEN];
      
for ii=1:length(A.SILCAT)
    if isnan(A.SILCAT(ii))==1 && isnan(A.PH_TOT(ii))==0 
        SILCAT_ca(ii)=silcat_nncanyonb_bit18(data(ii,:));
    else
        SILCAT_ca(ii)=A.SILCAT(ii);
    end
    if isnan(A.PHSPHT(ii))==1 && isnan(A.PH_TOT(ii))==0 
        PHSPHT_ca(ii)=phspht_nncanyonb_bit18(data(ii,:));
    else
        PHSPHT_ca(ii)=A.PHSPHT(ii);          
    end
end
   
SILCAT_ca=SILCAT_ca';
PHSPHT_ca=PHSPHT_ca';
clear data
      
%Caluclate missing TA values
for k=1:length(A.ALKALI)
    if isnan(A.ALKALI(k))==1 
       TA_ca(k)=A.SALNTY(k)*67;
    else
       TA_ca(k)=A.ALKALI(k);
    end
end
      
% Convert pH
[DATA,~,~]=CO2SYS(TA_ca(indph),A.PH_TOT(indph),1,3,A.SALNTY(indph),9,25,0,0,SILCAT_ca(indph),PHSPHT_ca(indph),1,10,1);
A.PH_TOT(indph)=DATA(:,37);
A.PH_TMP(:)=25;
clear DATA ii k PHSPHT_ca SILCAT_ca TA_ca indph

% Keep only flag=2 and flag=9 for TA and PH
A.PH_TOT(A.PH_TOT_FLAG_W~=2 & A.PH_TOT_FLAG_W~=6)=NaN;
A.PH_TOT_FLAG_W(A.PH_TOT_FLAG_W~=2 & A.PH_TOT_FLAG_W~=6)=9;
A.PH_TOT_FLAG_W(A.PH_TOT_FLAG_W==6)=2;
A.PH_TOT_FLAG_W(isnan(A.PH_TOT_FLAG_W))=9;
A.ALKALI(A.ALKALI_FLAG_W~=2 & A.ALKALI_FLAG_W~=6)=NaN;
A.ALKALI_FLAG_W(A.ALKALI_FLAG_W~=2 & A.ALKALI_FLAG_W~=6)=9;
A.ALKALI_FLAG_W(A.ALKALI_FLAG_W==6)=2;
A.ALKALI_FLAG_W(isnan(A.ALKALI_FLAG_W))=9;
A.TCARBN(A.TCARBN_FLAG_W~=2 & A.TCARBN_FLAG_W~=6)=NaN;
A.TCARBN_FLAG_W(A.TCARBN_FLAG_W~=2 & A.TCARBN_FLAG_W~=6)=9;
A.TCARBN_FLAG_W(A.TCARBN_FLAG_W==6)=2;
A.TCARBN_FLAG_W(isnan(A.TCARBN_FLAG_W))=9;

% Replace values from GLODAP with new values
ind=find(G2cruise==1046);
G2tdn(ind)=A.TDN;
G2phts25p0(ind)=A.PH_TOT;
G2phtsinsitutp(ind)=NaN;
G2talk(ind)=A.ALKALI;
G2tco2(ind)=A.TCARBN;
G2fco2(ind)=A.PCO2; % all NaN

G2tdnf(ind)=A.TDN_FLAG_W;
G2phts25p0f(ind)=A.PH_TOT_FLAG_W;
G2talkf(ind)=A.ALKALI_FLAG_W;
G2tco2f(ind)=A.TCARBN_FLAG_W;
G2fco2f(ind)=9;
G2phtsinsitutpf(ind)=9;
clearvars -except ind G2* expocode*

% Redo carbon calculations with new data
% Calculate "missing" SiO2 and PO4 using CANYONB
silicate=G2silicate;
phosphate=G2phosphate;
Si=CANYONB(G2year(isnan(silicate)),G2latitude(isnan(silicate)),G2longitude(isnan(silicate)),G2pressure(isnan(silicate)),G2temperature(isnan(silicate)),G2salinity(isnan(silicate)),G2oxygen(isnan(silicate)),'SiOH4');
silicate(isnan(silicate))=Si.SiOH4;
Phos=CANYONB(G2year(isnan(phosphate)),G2latitude(isnan(phosphate)),G2longitude(isnan(phosphate)),G2pressure(isnan(phosphate)),G2temperature(isnan(phosphate)),G2salinity(isnan(phosphate)),G2oxygen(isnan(phosphate)),'PO4');
phosphate(isnan(phosphate))=Phos.PO4;
clear Phos Si

% Calculate missing TA using S*67 approx.
talk=G2talk;
talk(isnan(talk))=G2salinity(isnan(talk)).*67;

% Find missing TCO2 values which can be calculated 
ind_misstco=find(isnan(G2tco2(ind))==1 & G2phts25p0f(ind)==2 & G2talkf(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
if length(ind_misstco)>=2*length(find(G2tco2f(ind)==2)) & length(find(G2tco2f(ind)==2))>0
    % If the conditions (described above) applies set all values to NaN
    % which later on can be replaced by calculated ones
    ind_misstco_new=find(G2phts25p0f(ind)==2 & G2talkf(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
    G2tco2f(ind(ind_misstco_new))=9;
    G2tco2(ind(ind_misstco_new))=NaN;
end

% Find missing TA values which can be calculated 
ind_misstalk=find(G2tco2f(ind)==2 & G2phts25p0f(ind)==2 & isnan(G2talk(ind))==1 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
if length(ind_misstalk)>=2*length(find(G2talkf(ind)==2)) & length(find(G2talkf(ind)==2))>0
    ind_misstalk_new=find(G2phts25p0f(ind)==2 & G2tco2f(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
    G2talkf(ind(ind_misstalk_new))=9;
    G2talk(ind(ind_misstalk_new))=NaN;
end
    
% Find missing PH_TOT values which can be calculated 
ind_misspHtot=find(G2tco2f(ind)==2 & isnan(G2phts25p0(ind))==1 & G2talkf(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
if length(ind_misspHtot)>=2*length(find(G2phts25p0f(ind)==2)) & length(find(G2phts25p0f(ind)==2))>0
    ind_misspHtot_new=find(G2tco2f(ind)==2 & G2talkf(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
    G2phts25p0f(ind(ind_misspHtot_new))=9;
    G2phts25p0(ind(ind_misspHtot_new))=NaN;
end
    
% Now use "final" indexes to replace missing carbon parameters - not
% considering pco2 for now
ind_misstco=ind(ind_misstco);
ind_misstalk=ind(ind_misstalk);
ind_misspHtot=ind(ind_misspHtot);

% Calculate missing parameters using CO2SYS
misstco=CO2SYS(G2talk(ind_misstco),G2phts25p0(ind_misstco),1,3,G2salinity(ind_misstco),25,G2temperature(ind_misstco),0,G2pressure(ind_misstco),silicate(ind_misstco),phosphate(ind_misstco),1,10,1);
G2tco2(ind_misstco)=misstco(:,2);
misstalk=CO2SYS(G2tco2(ind_misstalk),G2phts25p0(ind_misstalk),2,3,G2salinity(ind_misstalk),25,G2temperature(ind_misstalk),0,G2pressure(ind_misstalk),silicate(ind_misstalk),phosphate(ind_misstalk),1,10,1);
G2talk(ind_misstalk)=misstalk(:,1);
misspHtot=CO2SYS(G2talk(ind_misspHtot),G2tco2(ind_misspHtot),1,2,G2salinity(ind_misspHtot),G2temperature(ind_misspHtot),25,G2pressure(ind_misspHtot),0,silicate(ind_misspHtot),phosphate(ind_misspHtot),1,10,1);
G2phts25p0(ind_misspHtot)=misspHtot(:,18);

% Set flags for calculated values to 0 
G2tco2f(ind_misstco)=0;
G2talkf(ind_misstalk)=0;
flag0=find(~isnan(G2phts25p0(ind_misspHtot)));
G2phts25p0f(ind_misspHtot(flag0))=0;

% Missing pHinsitu at last as now calculated pHtot values can be used too
ind_misspHinsitu=find(isnan(G2phts25p0(ind))==0 & isnan(G2phtsinsitutp(ind))==1 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
ind_misspHinsitu=ind(ind_misspHinsitu);
misspHinsitu=CO2SYS(talk(ind_misspHinsitu),G2phts25p0(ind_misspHinsitu),1,3,G2salinity(ind_misspHinsitu),25,G2temperature(ind_misspHinsitu),0,G2pressure(ind_misspHinsitu),silicate(ind_misspHinsitu),phosphate(ind_misspHinsitu),1,10,1);
G2phtsinsitutp(ind_misspHinsitu)=misspHinsitu(:,18);

% Set pHinsitu flag to either 0 or 2 depending on pHtot
flag0=find(~isnan(G2phtsinsitutp(ind_misspHinsitu)) & G2phts25p0f(ind_misspHinsitu)~=2);
flag2=find(~isnan(G2phtsinsitutp(ind_misspHinsitu)) & G2phts25p0f(ind_misspHinsitu)==2);
flag9=find(G2phtsinsitutpf~=9 & isnan(G2phtsinsitutp)==1);
G2phtsinsitutpf(ind_misspHinsitu(flag2))=2;
G2phtsinsitutpf(ind_misspHinsitu(flag0))=0;
G2phtsinsitutpf(flag9)=9;

clearvars -except G2* talk silicate phosphate ind expocode*

% Also find missing pco2 values, i.e samples where at least two other carbon parameters are present; Always use DIC and TA if possible
ind_misspco=find(isnan(G2fco2(ind))==1 & G2talkf(ind)==2 & G2tco2f(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
ind_misspco_ta=find(isnan(G2fco2(ind))==1 & G2phts25p0f(ind)==2 & G2talkf(ind)==2 & G2tco2f(ind)~=2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);
ind_misspco_tco=find(isnan(G2fco2(ind))==1 & G2phts25p0f(ind)==2 & G2talkf(ind)~=2 & G2tco2f(ind)==2 & isnan(G2salinity(ind))==0 & isnan(G2oxygen(ind))==0);

ind_misspco=ind(ind_misspco);
ind_misspco_ta=ind(ind_misspco_ta);
ind_misspco_tco=ind(ind_misspco_tco);

% FCO2
misspco=CO2SYS(G2talk(ind_misspco),G2tco2(ind_misspco),1,2,G2salinity(ind_misspco),G2temperature(ind_misspco),20,G2pressure(ind_misspco),0,silicate(ind_misspco),phosphate(ind_misspco),1,10,1);
G2fco2(ind_misspco)=misspco(:,20);
misspco_ta=CO2SYS(G2talk(ind_misspco_ta),G2phts25p0(ind_misspco_ta),1,3,G2salinity(ind_misspco_ta),25,20,0,0,silicate(ind_misspco_ta),phosphate(ind_misspco_ta),1,10,1);
G2fco2(ind_misspco_ta)=misspco_ta(:,20);
misspco_tco=CO2SYS(G2phts25p0(ind_misspco_tco),G2tco2(ind_misspco_tco),3,2,G2salinity(ind_misspco_tco),25,20,0,0,silicate(ind_misspco_tco),phosphate(ind_misspco_tco),1,10,1);
G2fco2(ind_misspco_tco)=misspco_tco(:,20);

% Assign flag 0 to calculated data
G2fco2f(ind_misspco)=0;
G2fco2f(ind_misspco_ta)=0;
G2fco2f(ind_misspco_tco)=0;

% Set flag to 9 for missing pco2 values which could not be calculated due
% to missing silicate and phosphate, i.e. missing temperature
% etc...
na=find(isnan(G2fco2(ind)));
G2fco2f(ind(na))=9;

clear ind* miss* phosphate silicate talk


%% Update 33MW19910711 C14 
cd 'D:\GLODAPv2.2021\Updates\33MW19910711'
A=load('33MW19910711.mat');

% Sort the same way GLODAP data is sorted, i.e. first using station and
% then pressure
A=struct2table(A);
A=sortrows(A, [68 16]);

% Replace values
ind=find(G2cruise==336);
A.DELC14(A.DELC14_FLAG_W==3)=NaN;
A.C14ERR(A.DELC14_FLAG_W==3)=NaN;
A.DELC14_FLAG_W(A.DELC14_FLAG_W==3)=9;
A.DELC14_FLAG_W(isnan(A.DELC14))=9;
G2c14(ind)=A.DELC14;
G2c14err(ind)=A.C14ERR;
G2c14f(ind)=A.DELC14_FLAG_W;

clear ind A

%% Update 33RO20161119 C14 and C13 and get "real" btlnbr
cd 'D:\GLODAPv2.2021\Updates\33RO20161119'
A=load('33RO20161119.mat');

% Sort the same way GLODAP data is sorted, i.e. first using station and
% then pressure
A=struct2table(A);
A=sortrows(A, [68 16]);

% Replace values
ind=find(G2cruise==1045);
A.DELC14(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=NaN;
A.C14ERR(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=NaN;
A.DELC14_FLAG_W(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=9;
A.DELC14_FLAG_W(A.DELC14_FLAG_W==6)=2;
A.DELC13(A.DELC13_FLAG_W~=2 & A.DELC13_FLAG_W~=6)=NaN;
A.DELC13_FLAG_W(A.DELC13_FLAG_W~=2 & A.DELC13_FLAG_W~=6)=9;
A.DELC13_FLAG_W(A.DELC13_FLAG_W==6)=2;
G2c14(ind)=A.DELC14;
G2c14err(ind)=A.C14ERR;
G2c14f(ind)=A.DELC14_FLAG_W;
G2c13(ind)=A.DELC13;
G2c13f(ind)=A.DELC13_FLAG_W;
G2bottle(ind)=A.BTLNBR;

clear ind A

%% Correct bad CTDPRS of 33RO20131223 Station 5 castno 2

% Replace values --> ctdprs is "shifted" one bottle
ind=find(G2cruise==1042 & G2station==5 & G2cast==2);
holder=G2pressure(ind);
holderd=G2depth(ind);
for i=2:length(ind)
    G2pressure(ind(i))=holder(i-1);
    G2depth(ind(i))=holderd(i-1);
end

G2pressure(ind(1))=5804.0000;
G2depth(ind(1))=sw_dpth(5804.0000,G2latitude(ind(1)));
G2pressure(ind(3))=18.0000;
G2depth(ind(3))=sw_dpth(18.0000,G2latitude(ind(3)));

clear ind holder* i na


%% Save updated file

cd 'D:\GLODAPv2.2021\GLODAPv2.2020'
save('GLODAPv2.2021_Merged_Master_File_updated.mat')


