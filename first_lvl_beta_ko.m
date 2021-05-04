%% First level GLM model for spectral beta covariate on control subjects %%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b with SPM12 toolbox

% 1. Loads covariates, fmri scans, wake regressors, movement regressors
% 2. Builds SPM batch with the given parameters
% 3. Computes beta coeffiecients for every subject and saves

%-------------------------------------------------------------------------%
% GLM prvni urovne beta kovariaty na kontrolach%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Testovano na Matlab 2020b s toolboxem SPM12.


% 1. Nacte kovariaty, fmri scany, wake regresory, regresory pohybu
% 2. Prida vse do SPM batch 
% 3. Vypocita beta koeficienty glm modelu kazdeho subjekta a ulozi
%% Initialization, data loading, Inicializace, nacteni dat

close all; clear all; clc;

% Dir of the data, slozka s daty
mainPath = '/hydra-db/hydra_io/vypocty/piorecky/simult_sleep_data_MRcorrect/control/';
files = dir(mainPath);
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
dirFlags(1,5)=0;
subFolders = files(dirFlags);

% EEG covariate, EEG kovariata
covs_Path='/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/kovariaty/Occipital/Ko/';
covs_files = dir([covs_Path '*noinf_mean.mat']); % alfa, beta, theta, delta

% Wake correction, Korekce pocatecniho bdeni
WPath = '/hydradb/hydra_io/vypocty/piorecky/kovariaty/ko/W_regres/';
Wfiles = dir([WPath '*.mat']);

% For each subject, pro kazdeho subjekta
for i= 1:length(subFolders)   
subjectFile = [mainPath subFolders(i).name '/'];
subFile = dir(subjectFile);
dirFlag = [subFile.isdir] & ~strcmp({subFile.name},'.') & ~strcmp({subFile.name},'..') & ~strcmp({subFile.name},'*SPM*')& ~strcmp({subFile.name},'*lvl*') ;
subFolder = subFile(dirFlag);

cesta=[subjectFile subFolder(1).name '/'];

% Covariate, kovariata
posReg = ~cellfun('isempty',strfind({covs_files.name},subFolders(i).name));
rPathATh = [covs_Path covs_files(find(posReg == 1)).name];
load(rPathATh) %alfa, beta, delta, theta

swar_files = dir(fullfile(cesta,'*swarbfc*nii'));
addpath('/usr/local/spm12')
spm('defaults', 'fmri')
spm_jobman('initcfg')

for t=1:size(swar_files,1)
    swar_full_path(t,1)=cellstr(fullfile(cesta,swar_files(t).name));
end

% Motion correction file, korekce pohybu: 
motion_file = dir(fullfile(cesta,'rp*txt'));
a_motion = importdata(fullfile(cesta,motion_file.name));

%Wake regresor, regresor bdění subjekta:
posW =  ~cellfun('isempty',strfind({Wfiles.name},subFolders(i).name));
rPathW = [WPath Wfiles(find(posW == 1)).name];
load(rPathW) %k_delta k_alfa k_beta k_theta
       
pos_cov =  ~cellfun('isempty',strfind({covs_files.name},subFolders(i).name));
nazevSubj=covs_files(find(pos_cov == 1)).name(1:20);
first_name=['first_lvl_beta_',nazevSubj];
first_path=fullfile(['/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/'],first_name);

% fMRI settings spec., nastaveni fmri parametru:
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(first_path);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = swar_full_path;

% Match length of covariate and fmri
% srovnani delek fmri a eeg kovariaty

regress(1).name = 'anatnoninfmean';
regress(2).name = 'wake';
regress(1).val = k_beta';
regress(2).val = w1W';

if length(k_beta) > max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans))
    
      
    diff = length(k_beta)-max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans));
    val_reg = k_beta(1:end-diff);
    
     if isempty(isnan(val_reg)) == 0
        val_reg(isnan(val_reg)) = 0;
     end
    regress(1).val = val_reg;
    
end    
if length(w1W) > max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans))    
   
    diff = length(w1W)-max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans));
    valW = w1W(1:end-diff);
    
     if isempty(isnan(valW)) == 0
        valW(isnan(valW)) = 0;
     end
    
    regress(2).val = valW;
   
end    

if length(k_beta) < max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans)) 

    val_reg = zeros(1,max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans)));
    val_reg(1:length(k_beta)) = k_beta;

    if isempty(isnan(val_reg)) == 0
        val_reg(isnan(val_reg)) = 0;
    end
    
    regress(1).val = val_reg;
end
  
if length(w1W) < max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans))

    valW = zeros(1,max(size(matlabbatch{1}.spm.stats.fmri_spec.sess.scans)));
    valW(1:length(w1W)) = w1W;

    if isempty(isnan(valW)) == 0
        valW(isnan(valW)) = 0;
    end

    regress(2).val = valW'; %valA'    
    
    
    
end

% Build SPM batch, sestaveni SPM batch
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = regress;

matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(cesta,motion_file.name)}; %matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = subFolders(i).name;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1; % here change contrast [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

% Save Results, ulozeni vysledku: 
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = [subFolders(i).name i];
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec(1).thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{4}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec(1).mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.png = true;

spm_jobman('run',matlabbatch)

close all;
end

