% Read SBML model with COBRA Tool before executing FBA
% Execute 'initCobraToolbox(false)', and 
% 'modelOri = readCbModel('NAME.xml')' on the command window
% to make and reserve MAT files on the work space.
% Save this file on the current folder( 'save modelOri' on the command window).
clear
load modelOri.mat

exchgFnd = contains(modelOri.rxns,'EX_',"IgnoreCase",false);
rw_idx = find(exchgFnd == 1);
exchgTbl = table(rw_idx,modelOri.rxnNames(exchgFnd),...
    modelOri.rxns(exchgFnd),...
    modelOri.lb(exchgFnd),...
    modelOri.ub(exchgFnd),...
    'VariableNames', {'rxnRow','rxnNames','rxns','lb','ub'});
%% Start FBA
% set equation row
bio_syn = 7;  % Biomass synthesis equation
sbt_r = 180;  % substrate comsume equation, Glucose
sbt = -10;    % substrate comsumption value
co2_r = 946;
co2_cons = 0;
atpm_r = 645; % ATPM	ATP maintenance requirement
atpm = 3.15;  % ATP maintenace value mmol/gDCW/h
o2_r = 366; % o2 consumption equation % EX_o2_e	O2 exchange
o2_cons = -1000;     % maximum o2 consumption value mmol/gDCW/h
%% Setting parameters need for glpk function
rxns = modelOri.rxns;
rxnNames = modelOri.rxnNames;
c = zeros(size(modelOri.c,1),1); % vector of 0 other than objective function(biomass or target)
matrix = modelOri.S; % coefficient matrix for FBA
b = modelOri.b; % zeros(number of etabolites) vector
l = modelOri.lb; % lower bound (number of reactions) vector
u = modelOri.ub; % upper bound (number of reactions) vector
%% *Set Upper bound & Lower bound fluxes*
l(8,1)=0;
u(8,1)=0;     % 'BIOMASS_Ec_iJO1366_core_53p95M''E. coli biomass objective function (iJO1366) - core - with 53.95 GAM estimate'
l(101,1)= 0;
u(101,1)= 0;    % EX_fe3_e	Fe3+ exchange
l(1042,1)=0;
u(1042,1)=0;    % 'F6PA'	'Fructose 6-phosphate aldolase'
l(2554,1) = 0;
u(2554,1) = 0; % 'GLCtex_copy2'	'R_GLCtex_copy2' glc__D[e]  -> glc__D[p]
lb_grth = 0.05; % lowerbound of growth
l(bio_syn,1)= lb_grth;
l(o2_r,1) = o2_cons;
l(co2_r,1) = co2_cons; % 'CO2 exchange' 'EX_co2_e'
l(atpm_r,1)=atpm; % ATP maintenance requirement: atp[c] + h2o[c]  -> adp[c] + h[c] + pi[c]
%% Substrate exchange
l(sbt_r,1)=sbt;
u(sbt_r,1)=sbt;
%% csenseGLPK:
mtbnumbr = size(modelOri.metNames,1);
csense=repmat('E',mtbnumbr,1);
csenseGLPK = repmat('S',mtbnumbr,1);
%% Max O2_uptake set
c(bio_syn,1) = 1;
% Use only PTS
u(2552,1) = 0; % 'D-glucose transport in via proton symport (periplasm)'
Bmax = glpk(c,matrix, b, l, u, csenseGLPK); % Execute FBA for Biomass Max
l(o2_r,1) = Bmax(o2_r,1);
u(2552,1) = 1000;
%% Iterative objective function set
for i = 1:size(exchgTbl,1)
    trgt = exchgTbl.rxnRow(i,1);    % Target equation row number
    c(trgt,1) = 1;   % Maximize trgt
    c(bio_syn,1) = 10^-5;
    % In case of using only PTS
    u(2552,1) = 0; % 'D-glucose transport in via proton symport (periplasm)'
    TrgtMax_PTS = glpk(c,matrix, b, l, u, csenseGLPK); % Execute FBA for targetMax
    u(2552,1) = 1000; % 'D-glucose transport in via proton symport (periplasm)'
    % In case of using only Hxiokinase
    u(2551,1) = 0; % 'D-glucose transport via PEP:Pyr PTS (periplasm)'
    TrgtMax_HXK = glpk(c,matrix, b, l, u, csenseGLPK); % Execute FBA for targetMax
    u(2551,1) = 1000;
    if TrgtMax_PTS(trgt,1) >= TrgtMax_HXK(trgt,1)
        TrgtMax(:,i) = TrgtMax_PTS;
    else
        TrgtMax(:,i) = TrgtMax_HXK;
    end
    %     TrgtMax(:,i) = glpk(c,matrix, b, l, u, csenseGLPK); % Execute FBA for targetMax
    TrgtFlux(i,1) = TrgtMax(trgt,i);
    o2consFlux(i,1) = TrgtMax(o2_r,i);
    c = zeros(size(modelOri.c,1),1);
end
%% Data arrangement
exchgTbl.lb = l(exchgFnd);
exchgTbl.ub = u(exchgFnd);
TrgtFlux = table(TrgtFlux,'Variablenames',{'TrgtFlux'});
o2consFlux = table(o2consFlux,'Variablenames',{'o2consFlux'});
exchgTbl = [exchgTbl,TrgtFlux,o2consFlux];
toDelete = exchgTbl.TrgtFlux < 1.0; 
exchgTbl(toDelete,:) = [];
TrgtFlux(toDelete,:) = [];
TrgtMax(:,toDelete) = [];
toDelete = exchgTbl.lb < 0; % remove as a target such as h, h2o, fe2, fe3 etc.
exchgTbl(toDelete,:) = [];
TrgtFlux(toDelete,:) =[];
TrgtMax(:,toDelete) =[];
% save info1
%% Prepare AERITH
% load info1.mat
% select Exchange target
% As an example, succinate production is performed.
% exchgTbl(8): 'Succinate exchange','EX_succ_e'

% for j = 1:size(exchgTbl,1)
j = 8;
exchgTrgt = j;
trgt = exchgTbl.rxnRow(exchgTrgt,1);    % Target equation row number
c = zeros(size(modelOri.c,1),1);
c(bio_syn,1) = 1;   % Maximize biomass
c(trgt,1) = 10^-5; % use to maximize growth and then a target(avoid indeterminate answers)
u(exchgTbl.rxnRow,1)=0;
u(trgt,1)=1000;
l(o2_r,1) = exchgTbl.o2consFlux(exchgTrgt,1);
u(o2_r,1) = exchgTbl.o2consFlux(exchgTrgt,1);
u(43,1)=1000;     % 'Succinate exchange'	'EX_succ_e'
u(95,1)=1000;     % 'EX_etoh_e'
u(119,1)=1000;     % 'Urea exchange'	'EX_uea_e'
u(125,1)=1000;  % 'EX_for_e'
u(167,1)=1000;  % '5-Methylthio-D-ribose exchange'	'EX_5mtr_e' Growth essential
u(223,1)=1000;    % EX_ac_e	Acetate exchange
u(227,1)=1000;    % EX_h2s_e	Hydrogen sulfide exchange
u(266,1)=1000;  % 'EX_lac__D_e'
u(290,1)=1000;  % EX_mal__L_e	L-Malate exchange
u(825,1)=1000;  % 'Citrate exchange'	'EX_cit_e'
u(946,1)=1000;  % 'CO2 exchange'	'EX_co2_e'
normalFBAresult = glpk(c,matrix, b, l, u, csenseGLPK);
%% Remove exchange reactions
posKO_reac = true(size(rxns,1),1);
posKO_reac = posKO_reac-exchgFnd;
%% Remove other transport reactions
for i = 1:size(rxns,1)
    if posKO_reac(i,1) == 1
        if (isempty(strfind(rxns{i},'tipp')) == 0 && strfind(rxns{i},'tipp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'tex')) == 0 && strfind(rxns{i},'tex')+2 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'tpp')) == 0 && strfind(rxns{i},'tpp')+2 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t2pp')) == 0 && strfind(rxns{i},'t2pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t2rpp')) == 0 && strfind(rxns{i},'t2rpp')+4 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t3ipp')) == 0 && strfind(rxns{i},'t3ipp')+4 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t7pp')) == 0 && strfind(rxns{i},'t7pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t3pp')) == 0 && strfind(rxns{i},'t3pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t8pp')) == 0 && strfind(rxns{i},'t8pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t2_2pp')) == 0 && strfind(rxns{i},'t2_2pp')+5 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t4pp')) == 0 && strfind(rxns{i},'t4pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t2_3pp')) == 0 && strfind(rxns{i},'t2_3pp')+5 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t2ipp')) == 0 && strfind(rxns{i},'t2ipp')+4 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t6pp')) == 0 && strfind(rxns{i},'t6pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t3')) == 0 && strfind(rxns{i},'t3')+1 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'trpp')) == 0 && strfind(rxns{i},'trpp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t6_2pp')) == 0 && strfind(rxns{i},'t6_2pp')+5 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'texi')) == 0 && strfind(rxns{i},'texi')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'t1pp')) == 0 && strfind(rxns{i},'t1pp')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'tppi')) == 0 && strfind(rxns{i},'tppi')+3 == size(rxns{i},2)) ...
                || (isempty(strfind(rxns{i},'tiex')) == 0 && strfind(rxns{i},'tiex')+3 == size(rxns{i},2))
            posKO_reac(i,1) = 0;
        end
    end
end
bffr = 0.90; % r value for calculation flexibility
normalFBAresult(abs(normalFBAresult)<10^-5,1) = 0;
TrgtMax(abs(TrgtMax(:,exchgTrgt))< 10^-5,exchgTrgt) = 0;
rct_row = (1:size(rxns,1))';
AnalysisTble = table(rct_row,rxns,rxnNames,posKO_reac);
chg_rate = (abs(TrgtMax(:,exchgTrgt))-abs(normalFBAresult))./abs(normalFBAresult);
chg_rate(isnan(chg_rate))=0; % normalFBAresult: 0, TrgtMax: 0
chg_rate(isinf(chg_rate))=10000; % normalFBAresult: 0
chg_rate(chg_rate<min(chg_rate)*bffr,1) = min(chg_rate);
absFBArslt = abs(normalFBAresult);
AnalysisTble = [AnalysisTble table(chg_rate) table(absFBArslt) table(normalFBAresult)...
    table(TrgtMax(:,exchgTrgt),'VariableNames',{'TrgtMax'})];
AnalysisTble = sortrows(AnalysisTble,{'chg_rate','absFBArslt'},{'ascend','descend'});
AnalysisTble = sortrows(AnalysisTble,'rct_row');

AnalysisTble.posKO_reac(1808,1) = 0; % 'NAt3_1p5pp'	'Sodium proton antiporter (H:NA is 1.5) (periplasm)'
AnalysisTble.posKO_reac(2553,1) = 0; % GLCtex_copy1	R_GLCtex_copy1
AnalysisTble.posKO_reac(2554,1) = 0; % GLCtex_copy2	R_GLCtex_copy2
AnalysisTble = sortrows(AnalysisTble,{'chg_rate','absFBArslt'},{'ascend','descend'});
%% AERITH: Algorithm of Efficiently Reaction Identifications for a Target compound with High productivity
maxKO = 100; % maximum numbers of knock out reactions
l_new = l;
u_new = u;
l_sko = l_new;
u_sko = u_new;
FBArslt_KO = zeros(size(rxns,1),maxKO);
TrgtMaxFBArslt_KO = zeros(size(rxns,1),maxKO);
KOgrwth = zeros(maxKO,1);
KOtrgt = zeros(maxKO,1);
optKOset = [];
KOset = [];
for k = 1 : maxKO
    for i = 1:size(rxns,1)
        %         if AnalysisTble.posKO_reac(i,1) == 1 && AnalysisTble.absFBArslt(i,1) < 1000
        if AnalysisTble.posKO_reac(i,1) == 1 && AnalysisTble.absFBArslt(i,1) < 100
            l_sko(AnalysisTble.rct_row(i,1),1) = 0;
            u_sko(AnalysisTble.rct_row(i,1),1) = 0;
            c = zeros(size(modelOri.c,1),1); % vector of 0 other than objective function(biomass or target)
            c(bio_syn,1) = 1;
            c(trgt,1) = 10^-5;
            Grth_chk = glpk(c,matrix, b, l_sko, u_sko, csenseGLPK); % biomass->max
            AnalysisTble.posKO_reac(i,1) = 0; 
            if Grth_chk(bio_syn,1) >= lb_grth  % Checking if possible to grow even if knockout
                break
            else
                l_sko = l_new;
                u_sko = u_new;
            end
        end
    end
    KOset = vertcat(KOset, AnalysisTble(i,1:3));
    FBArslt_KO(:,k) = Grth_chk;
    FBArslt_KO(abs(FBArslt_KO(:,k)) < 10^-5,k) = 0;
    absFBArslt = abs(FBArslt_KO(:,k));
    KOgrwth(k,1) = FBArslt_KO(bio_syn,k);
    
    l_grth = l_sko;
    c = zeros(size(modelOri.c,1),1); % vector of 0 other than objective function(biomass or target)
    c(trgt,1) = 1;   % Maximize trgt
    c(bio_syn,1) = 10^-5;
    TrgtMaxFBArslt_KO(:,k) = glpk(c,matrix, b, l_grth, u_sko, csenseGLPK); % target ->max
    KOtrgt(k,1) = FBArslt_KO(trgt,k);
    chg_rate = (abs(TrgtMaxFBArslt_KO(:,k))-abs(FBArslt_KO(:,k)))./abs(FBArslt_KO(:,k));
    chg_rate(isnan(chg_rate))=0;
    chg_rate(isinf(chg_rate))=10000;
    chg_rate(chg_rate<min(chg_rate)*bffr,1) = min(chg_rate);
    AnalysisTble = AnalysisTble(:,1:4);
    AnalysisTble = sortrows(AnalysisTble,"rct_row");
    AnalysisTble = [AnalysisTble table(chg_rate) table(absFBArslt)];
    AnalysisTble = sortrows(AnalysisTble,{'chg_rate','absFBArslt'},{'ascend','descend'});
    l_new = l_sko;
    u_new = u_sko;
end
optGrthTrgt = [(FBArslt_KO(bio_syn,:))' (FBArslt_KO(trgt,:))'];
optKOrslt = [KOset table(optGrthTrgt)];

exchg_rslt = FBArslt_KO(exchgTbl.rxnRow,:);
exchg_rslt = table(exchg_rslt,'VariableNames',{'Rslt_Flux'});
exchg_rslt = [exchgTbl(:,1:3),exchg_rslt];
rootname = 'result';
extension = '.mat';
filename = [rootname, num2str(exchgTrgt), extension];
save(filename)