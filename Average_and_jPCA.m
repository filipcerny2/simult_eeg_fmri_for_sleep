%% Averaging and jPCA EEG dimension reduction method %%
        % For anatomically reduced covariates 
        
% Author:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b.

% 1. loads spectral EEG data of anatomical region
% 2. Reduces spectrum to 1D vector by averaging power across EEG channels.
% 3. Reduced spectrum to 1D vector by group PCA.
%  3.1 Subtracts mean value of every subject.
%  3.2 Lines data of each subject into one matrix.
%  3.3 PCA is performed and the first component is isolated.
%  3.4 Data is split to individual subjects and saved. 

%-------------------------------------------------------------------------%
% Metoda prumerovani a jPCA na antomicky redukovaných kovariátach %

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Testováno na Matlab 2020b.

% 1. Skript nacte vytvorena EEG spektra anatomickych regionu.
% 2. Spektra jsou redukovana na 1D vektor prumerovanim vykonu napric kanaly
% v regionu.
% 3. Spektra jsou redukovana na 1D vektor metodou spolecne PCA 
%   3.1. Je odectena stredni hodnota spektra kazdeho subjekta
%   3.2. Spektra subjektu jsou poskladany za sebe
%   3.3. Je provedena PCA a izolovana prvni komponenta (1. sloupec)
%   3.4. Data jsou opet rozdelena na subjekty a ulozena

%% Initialization, Inicializace a nacteni datoveho souboru
clc; clear all; close all;

gen_Path='D:\Marek_P\Filip\Anat\'; % nazvy slozek na anatomicke regiony
nazvy_slozek=dir(gen_Path);
folders=cell(6:1);
for k=3:8
    folders{k-2}=nazvy_slozek(k).name;
end

FilePath_ko=('D:\Marek_P\Filip\Anat\Frontal\Ko\');
FileName_ko = dir([FilePath_ko '*.mat']);
nFiles_ko = length(FileName_ko);

FilePath_pac=('D:\Marek_P\Filip\Anat\Frontal\Pac\');
FileName_pac = dir([FilePath_pac '*.mat']);
nFiles_pac = length(FileName_pac);


%% Averaging method, Prumerovani vykonu 

for e=1:length(folders)
    for i=1:nFiles_ko % kontroly
        load(['D:\Marek_P\Filip\Anat\' folders{e} '\Ko\' FileName_ko(i).name]);
        
        pwr.alfa=mean(pwr.alfa.powspctrm,2); % prumerovani
        pwr.beta=mean(pwr.beta.powspctrm,2);
        pwr.theta=mean(pwr.theta.powspctrm,2);
        pwr.delta=mean(pwr.delta.powspctrm,2);
        save(['D:\Marek_P\Filip\Anat_mean\' folders{e} '\Ko\' FileName_ko(i).name],'pwr','-v7.3');
    end

    for i=1:nFiles_pac % pacienti
        load(['D:\Marek_P\Filip\Anat\' folders{e} '\Pac\' FileName_pac(i).name]);
        
        pwr.alfa=mean(pwr.alfa.powspctrm,2); % prumerovani
        pwr.beta=mean(pwr.beta.powspctrm,2);
        pwr.theta=mean(pwr.theta.powspctrm,2);
        pwr.delta=mean(pwr.delta.powspctrm,2);
        save(['D:\Marek_P\Filip\Anat_mean\' folders{e} '\Pac\' FileName_pac(i).name],'pwr','-v7.3');
    end
end

%% Group PCA, jPCA spektralniho vykonu

mat_alfa=cell(1,8); 
mat_beta=cell(1,8);
mat_delta=cell(1,8);
mat_theta=cell(1,8);

for e=1:length(folders)
    for i=1:nFiles_ko % kontroly
        load(['D:\Marek_P\Filip\Anat\' folders{e} '\Ko\' FileName_ko(i).name]);
        
        % Subtract mean value, odecist stredni hodnotu
        
        mat_alfa1 = pwr.alfa.powspctrm;
        mat_beta1 = pwr.beta.powspctrm;
        mat_delta1 = pwr.delta.powspctrm;
        mat_theta1 = pwr.theta.powspctrm;
    
        str_alfa = std(mat_alfa1); % stredni hodnoty spekter
        str_beta = std(mat_beta1);
        str_delta = std(mat_delta1);
        str_theta = std(mat_theta1);
        
        for j= 1:size(mat_alfa1,1)  % odecteni strednich hodnot
            mat_alfa2(j,:) = mat_alfa1(j,:)-str_alfa;
            mat_beta2(j,:) = mat_beta1(j,:)-str_beta;
            mat_delta2(j,:) = mat_delta1(j,:)-str_delta;
            mat_theta2(j,:) = mat_theta1(j,:)-str_theta;
        end
        
        mat_alfa{1,i}=num2cell(mat_alfa2);
        mat_beta{1,i}=num2cell(mat_beta2);
        mat_delta{1,i}=num2cell(mat_delta2);
        mat_theta{1,i}=num2cell(mat_theta2);
       
        clear mat_alfa1 mat_beta1 mat_delta1 mat_theta1 mat_alfa2 mat_beta2 mat_delta2 mat_theta2;
        
    end
    
    % Ordering of the subjects in one matrix, Poskladani subjektu za sebe
    
    delky_ko=[length(mat_alfa{1});length(mat_alfa{2});length(mat_alfa{3});length(mat_alfa{4});...
        length(mat_alfa{5});length(mat_alfa{6});length(mat_alfa{7});length(mat_alfa{8})];
    comp_mat_alfa=cell2mat([mat_alfa{1};mat_alfa{2};mat_alfa{3};mat_alfa{4};...
            mat_alfa{5};mat_alfa{6};mat_alfa{7};mat_alfa{8}])';
    comp_mat_beta=cell2mat([mat_beta{1};mat_beta{2};mat_beta{3};mat_beta{4};...
            mat_beta{5};mat_beta{6};mat_beta{7};mat_beta{8}])';
    comp_mat_delta=cell2mat([mat_delta{1};mat_delta{2};mat_delta{3};mat_delta{4};...
            mat_delta{5};mat_delta{6};mat_delta{7};mat_delta{8}])';
    comp_mat_theta=cell2mat([mat_theta{1};mat_theta{2};mat_theta{3};mat_theta{4};...
            mat_theta{5};mat_theta{6};mat_theta{7};mat_theta{8}])';
     
    % PCA analysis, Provedeni PCA analyzy a izolace prvni komponenty   
         
    mat_alfa_pca=pca(comp_mat_alfa);
    jpca.alfa=mat_alfa_pca(:,1);
    
    mat_beta_pca=pca(comp_mat_beta);
    jpca.beta=mat_beta_pca(:,1);
    
    mat_delta_pca=pca(comp_mat_delta);
    jpca.delta=mat_delta_pca(:,1);
    
    mat_theta_pca=pca(comp_mat_theta);
    jpca.theta=mat_theta_pca(:,1);

    % Division back to subjects, rozdeleni zpet na jednotlive subjekty
    
    m=0;
    for r=1:nFiles_ko
        pca_cov.alfa=jpca.alfa((1+m):(delky_ko(r)+m));
        pca_cov.beta=jpca.beta((1+m):(delky_ko(r)+m));
        pca_cov.delta=jpca.delta((1+m):(delky_ko(r)+m));
        pca_cov.theta=jpca.theta((1+m):(delky_ko(r)+m));
        m=sum(delky_ko(1:r));
        save(['D:\Marek_P\Filip\Anat_pca\' folders{e} '\Ko\' FileName_ko(r).name],'pca_cov','-v7.3');
    end
    
    clear pca_cov jpca pwr comp_mat_alfa comp_mat_beta comp_mat_delta comp_mat_theta;
    
    
    % Patients, Pacienti
    
    for i=1:nFiles_pac
        load(['D:\Marek_P\Filip\Anat\' folders{e} '\Pac\' FileName_pac(i).name]);
        
        % Remove mean, odecist stredni hodnotu
        
        mat_alfa1 = pwr.alfa.powspctrm;
        mat_beta1 = pwr.beta.powspctrm;
        mat_delta1 = pwr.delta.powspctrm;
        mat_theta1 = pwr.theta.powspctrm;
    
        str_alfa = std(mat_alfa1); 
        str_beta = std(mat_beta1);
        str_delta = std(mat_delta1);
        str_theta = std(mat_theta1);
        
        for j= 1:size(mat_alfa1,1)
            mat_alfa2(j,:) = mat_alfa1(j,:)-str_alfa;
            mat_beta2(j,:) = mat_beta1(j,:)-str_beta;
            mat_delta2(j,:) = mat_delta1(j,:)-str_delta;
            mat_theta2(j,:) = mat_theta1(j,:)-str_theta;
        end
        
        mat_alfa{1,i}=num2cell(mat_alfa2);
        mat_beta{1,i}=num2cell(mat_beta2);
        mat_delta{1,i}=num2cell(mat_delta2);
        mat_theta{1,i}=num2cell(mat_theta2);
       
        clear mat_alfa1 mat_beta1 mat_delta1 mat_theta1 mat_alfa2 mat_beta2 mat_delta2 mat_theta2;
        
    end
    
    % Ordering of the subjects in one matrix, Poskladani subjektu za sebe
    
    delky_pac=[length(mat_alfa{1});length(mat_alfa{2});length(mat_alfa{3});length(mat_alfa{4});...
        length(mat_alfa{5});length(mat_alfa{6});length(mat_alfa{7});length(mat_alfa{8})];
    comp_mat_alfa=cell2mat([mat_alfa{1};mat_alfa{2};mat_alfa{3};mat_alfa{4};...
            mat_alfa{5};mat_alfa{6};mat_alfa{7};mat_alfa{8}])';
    comp_mat_beta=cell2mat([mat_beta{1};mat_beta{2};mat_beta{3};mat_beta{4};...
            mat_beta{5};mat_beta{6};mat_beta{7};mat_beta{8}])';
    comp_mat_delta=cell2mat([mat_delta{1};mat_delta{2};mat_delta{3};mat_delta{4};...
            mat_delta{5};mat_delta{6};mat_delta{7};mat_delta{8}])';
    comp_mat_theta=cell2mat([mat_theta{1};mat_theta{2};mat_theta{3};mat_theta{4};...
            mat_theta{5};mat_theta{6};mat_theta{7};mat_theta{8}])';
     
    % PCA analysis, Provedeni PCA analyzy a izolace prvni komponenty     
         
    mat_alfa_pca=pca(comp_mat_alfa);
    jpca.alfa=mat_alfa_pca(:,1);
    
    mat_beta_pca=pca(comp_mat_beta);
    jpca.beta=mat_beta_pca(:,1);
    
    mat_delta_pca=pca(comp_mat_delta);
    jpca.delta=mat_delta_pca(:,1);
    
    mat_theta_pca=pca(comp_mat_theta);
    jpca.theta=mat_theta_pca(:,1); 
    
    % Division back to subjects, rozdeleni zpet na jednotlive subjekty
    
    m=0;
    for r=1:nFiles_pac
        pca_cov.alfa=jpca.alfa((1+m):(delky_pac(r)+m));
        pca_cov.beta=jpca.beta((1+m):(delky_pac(r)+m));
        pca_cov.delta=jpca.delta((1+m):(delky_pac(r)+m));
        pca_cov.theta=jpca.theta((1+m):(delky_pac(r)+m));
        m=sum(delky_pac(1:r));
        save(['D:\Marek_P\Filip\Anat_pca\' folders{e} '\Pac\' FileName_pac(r).name],'pca_cov','-v7.3');
    end 
end