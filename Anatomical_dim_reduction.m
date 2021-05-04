
%% Anatomical dimension reduction of hdEEG space and spectral analysis%%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b with Fieldtrip toolbox

% 1. Loads EEG data in Fieldtrip toolbox format.
% 2. Data is split into regions with EGI_electrode_map.m script  
% 3. Spectral anylsis is performed into 4 spectral bands. 
% 4. Bands are normalize for inter subject comparison. 
%-------------------------------------------------------------------------%
% Anatomicka redukce dimenze hdEEG prostoru a spektralni analyza%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b.

% 1. Skript nacte EEG data pacientu a kontrol.
% 2. Data jsou rozdelena na anatomicke regiony dle elektrod helmy EGI256.  
% 3. Je provedena frekvencni analyza dat na 4 pasma. 
% 4. Pasma jsou normalizovany pro moznost porovnani.

%% Initialization and data loading, Inicializace a nacteni toolboxu a dat

clc; clear all; close all;

addpath D:\Vlasta\GIT\fieldtrip\ % fieldtrip toolbox
ft_defaults; % inicialization of fieldtrip

% load electrode_map and subject shifts
load tr_kont.mat; load tr_pac.mat; load EGI_electrode_map.m;

tr_kont= orderfields(tr_kont); % order shifts
tr_pac = orderfields(tr_pac); 

gen_Path='D:\Marek_P\Filip\Anat\'; % folders for anatomical regions
nazvy_slozek=dir(gen_Path);
folders=cell(6:1);
for k=3:8
    folders{k-2}=nazvy_slozek(k).name;
end

FilePath_ko=('D:\Marek_P\Filip\sekundy\Ko\'); % load control subjects,  nacteni dat kontrol
FileName_ko = dir([FilePath_ko '*.mat']);
nFiles_ko = length(FileName_ko);

FilePath_pac=('D:\Marek_P\Filip\sekundy\Pac\'); % load patients, nacteni dat pacientu
FileName_pac = dir([FilePath_pac '*.mat']);
nFiles_pac = length(FileName_pac);

ko=fieldnames(tr_kont); 
pac=fieldnames(tr_pac);


%% Anatomical dimension reduction and spectral analysis 
% Anatomicka redukce dimenze s frekvencni analyzou

for j=1:nFiles_ko                               % constrols, pro skupinu kontrol
    
    load([FilePath_ko FileName_ko(j).name]);    
    posun=round(tr_kont.(ko{j})/250);           
    vekt=zeros(1,length(data_sek1.trial));      
    vekt(1,posun:end)=1;
    cfg=[];
    cfg.trials=logical(vekt);
    data_sek1=ft_selectdata(cfg,data_sek1);     % cut to fmri size, zkrati na delku fmri
     
    % Region dimension reduction, Redukovani dimenze EEG prostoru na vybrany region
    
    % frontalni oblast
    cfg=[];
    cfg.channel=(el_map.frontal)';
    anat_ko.frontal=ft_preprocessing(cfg,data_sek1);
        
    % parietalni oblast
    cfg=[];
    cfg.channel=(el_map.parietal)';
    anat_ko.parietal=ft_preprocessing(cfg,data_sek1);
    
    % okcipitalni oblast
    cfg=[];
    cfg.channel=(el_map.occip)';
    anat_ko.occipital=ft_preprocessing(cfg,data_sek1);
    
    % temporalni oblast
    cfg=[];
    cfg.channel=(el_map.temporal)';
    anat_ko.temporal=ft_preprocessing(cfg,data_sek1);
    
    % leva hemisfera
    cfg=[];
    cfg.channel=(el_map.left_hem)';
    anat_ko.left_hem=ft_preprocessing(cfg,data_sek1);
    
    % prava hemisfera
    cfg=[];
    cfg.channel=(el_map.right_hem)';
    anat_ko.right_hem=ft_preprocessing(cfg,data_sek1);
    anat_ko=orderfields(anat_ko);
    NAnat=fieldnames(anat_ko); % pocet regionu
    
    % freq analysis, frekvencni analyza na alfa, beta, delta theta v regionech
    for i=1:length(NAnat)
    
        cfg = [];                
        cfg.trials       = 'all';     % parametry frekvencni analyzy
        cfg.output       = 'pow';
        cfg.method       = 'mtmfft';
        cfg.taper        = 'dpss';
        cfg.keeptrials   = 'yes';
    
        % DELTA
        cfg.foi          = 2;
        cfg.tapsmofrq    = 2;
        delta = ft_freqanalysis(cfg,anat_ko.(NAnat{i}));
    
        % THETA
        cfg.foi          = 6;
        cfg.tapsmofrq    = 2;
        theta = ft_freqanalysis(cfg,anat_ko.(NAnat{i}));
    
        % ALFA
        cfg.foi          = 10.5;
        cfg.tapsmofrq    = 2.5;
        alfa = ft_freqanalysis(cfg,anat_ko.(NAnat{i}));
        alfa.freq=10.5;
    
        % BETA
        cfg.foi          = 16.5;
        cfg.tapsmofrq    = 3.5;
        beta = ft_freqanalysis(cfg,anat_ko.(NAnat{i}));
        beta.freq=16.5;
    
        suma = sum(delta.powspctrm,1) ... % normalizace spekter
             + sum(theta.powspctrm,1) ...
             + sum(alfa.powspctrm,1) ...
             + sum(beta.powspctrm,1);
    
        GlobalPower = suma;
        suma = zeros(1,3);
         
         deltaPA = bsxfun(@rdivide,delta.powspctrm,GlobalPower); % normalization, normalizace
         thetaPA = bsxfun(@rdivide,theta.powspctrm,GlobalPower);
         alfaPA = bsxfun(@rdivide,alfa.powspctrm,GlobalPower);
         betaPA = bsxfun(@rdivide,beta.powspctrm,GlobalPower);
         
         delta.powspctrm  = deltaPA;
         theta.powspctrm  = thetaPA;
         alfa.powspctrm  = alfaPA;
         beta.powspctrm = betaPA;    
         
         pwr.delta=delta;
         pwr.theta=theta;
         pwr.alfa=alfa;
         pwr.beta=beta;
         
         Filename_new=sprintf('ko%d',j);
         save(['D:\Marek_P\Filip\Anat\' folders{1,i} '\Ko\'  Filename_new],'pwr','-v7.3');
         clear pwr alfa beta delta theta GlobalPower deltaPA alfaPA betaPA thetaPA;
    
    end   
end

% pacienti

for j=1:nFiles_pac                              % patients, pro skupinu pacientu
    
    load([FilePath_pac FileName_pac(j).name]);  
    posun=round(tr_pac.(pac{j})/250);           
    vekt=zeros(1,length(data_sek2.trial));      
    vekt(1,posun:end)=1;
    cfg=[];
    cfg.trials=logical(vekt);
    data_sek2=ft_selectdata(cfg,data_sek2);     % cuts to fmri size, zkrati na delku fmri
     
    % Redukovani dimenze EEG prostoru na vybrany region
    
    % frontalni oblast
    cfg=[];
    cfg.channel=(el_map.frontal)';
    anat_pac.frontal=ft_preprocessing(cfg,data_sek2);
        
    % parietalni oblast
    cfg=[];
    cfg.channel=(el_map.parietal)';
    anat_pac.parietal=ft_preprocessing(cfg,data_sek2);
    
    % okcipitalni oblast
    cfg=[];
    cfg.channel=(el_map.occip)';
    anat_pac.occipital=ft_preprocessing(cfg,data_sek2);
    
    % temporalni oblast
    cfg=[];
    cfg.channel=(el_map.temporal)';
    anat_pac.temporal=ft_preprocessing(cfg,data_sek2);
    
    % leva hemisfera
    cfg=[];
    cfg.channel=(el_map.left_hem)';
    anat_pac.left_hem=ft_preprocessing(cfg,data_sek2);
    
    % prava hemisfera
    cfg=[];
    cfg.channel=(el_map.right_hem)';
    anat_pac.right_hem=ft_preprocessing(cfg,data_sek2);
    anat_pac=orderfields(anat_pac);
    NAnat=fieldnames(anat_pac);
    
    % Frekvencni analyza na spektra alfa beta delta theta v regionech
    for i=1:length(NAnat)
    
        cfg = [];                
        cfg.trials       = 'all';     %viz help (nevybira dle podminky)
        cfg.output       = 'pow';
        cfg.method       = 'mtmfft';
        cfg.taper        = 'dpss';
        cfg.keeptrials   = 'yes';
    
        % DELTA
        cfg.foi          = 2;
        cfg.tapsmofrq    = 2;
        delta = ft_freqanalysis(cfg,anat_pac.(NAnat{i}));
    
        % THETA
        cfg.foi          = 6;
        cfg.tapsmofrq    = 2;
        theta = ft_freqanalysis(cfg,anat_pac.(NAnat{i}));
    
        % ALFA
        cfg.foi          = 10.5;
        cfg.tapsmofrq    = 2.5;
        alfa = ft_freqanalysis(cfg,anat_pac.(NAnat{i}));
        alfa.freq=10.5;
    
        % BETA
        cfg.foi          = 16.5;
        cfg.tapsmofrq    = 3.5;
        beta = ft_freqanalysis(cfg,anat_pac.(NAnat{i}));
        beta.freq=16.5;
    
        suma = sum(delta.powspctrm,1) ...
             + sum(theta.powspctrm,1) ...
             + sum(alfa.powspctrm,1) ...
             + sum(beta.powspctrm,1);
    
        GlobalPower = suma;
        suma = zeros(1,3);
         
         deltaPA = bsxfun(@rdivide,delta.powspctrm,GlobalPower);
         thetaPA = bsxfun(@rdivide,theta.powspctrm,GlobalPower);
         alfaPA = bsxfun(@rdivide,alfa.powspctrm,GlobalPower);
         betaPA = bsxfun(@rdivide,beta.powspctrm,GlobalPower);
         
         delta.powspctrm  = deltaPA;
         theta.powspctrm  = thetaPA;
         alfa.powspctrm  = alfaPA;
         beta.powspctrm = betaPA;    
         
         pwr.delta=delta;
         pwr.theta=theta;
         pwr.alfa=alfa;
         pwr.beta=beta;
         
         Filename_new=sprintf('pac%d',j);
         save(['D:\Marek_P\Filip\Anat\' folders{1,i} '\Pac\'  Filename_new],'pwr','-v7.3');
         clear pwr alfa beta delta theta GlobalPower deltaPA alfaPA betaPA thetaPA;
    
    end   
end


