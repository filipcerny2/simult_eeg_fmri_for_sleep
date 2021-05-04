%% Second level GLM model for spectral beta covariate on control subjects %%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b with SPM12 toolbox

% 1. Loads computed beta values from first level of every subject
% 2. Statistically computes similarly activated regions across subjects
% 3. Maps significantly activated regions on a common statistical
%  parametrical map.


%-------------------------------------------------------------------------%
% GLM druhe urovne beta kovariaty na kontrolach%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Testovano na Matlab 2020b s toolboxem SPM12.


% 1. Nacte vypocitane beta hodnoty z prvni urovne statistiky
% 2. Statisticky urci stejne aktivovane regiony napric subjekty. 
% 3. Signifikantni regiony jsou zapsany na statisticko parametrickou mapu.
%% Initialization, data loading, Inicializace, nacteni dat

clc; clear all; close all; 

addpath('/usr/local/spm12')
spm('defaults', 'fmri')
spm_jobman('initcfg')

% load dtaa from 1st level, nacteni dat z prvni urovne
first_res={ ...    
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02409_20190514/con_0001.nii'
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02411_20190507/con_0001.nii'
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02423_20181016/con_0001.nii'
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02426_20190205/con_0001.nii'
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02428_20190319/con_0001.nii'

'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02425_20181127/con_0001.nii'
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02427_20190305/con_0001.nii'
'/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/first_lvl/Occipital/beta/Ko/first_lvl_beta_PARS_P02430_20190411/con_0001.nii'
};

matlabbatch{1}.spm.stats.factorial_design.dir = cellstr('/hydra-db/hydra_io/vypocty/cerny/Anat_noninf_mean/second_lvl/beta/all_Ko/');
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = first_res;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Occ_beta_ko';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = 'Occ_beta_ko';
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = 'none';
matlabbatch{4}.spm.stats.results.conspec(1).thresh = 0.001;
matlabbatch{4}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{4}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec(1).mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.png = true;

spm('defaults', 'fmri')
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)


