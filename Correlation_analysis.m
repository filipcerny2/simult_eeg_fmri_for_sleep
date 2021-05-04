
%% Correlation analysis of anatomically reduced covariates for Mean and jPCA method %%

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b.

% 1. Loads anatomically reduced covariated in 1D vector format of the selected subject.
% 2. Vectors are ordered into one matrix. 
% 3. corcoeff function is used for obtaining the Pearson correalation coefficient matrix.

%-------------------------------------------------------------------------%

% Korelacni analyza anatomicky redukovanych kovariat pro metodu redukce dimenze jPCA %

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Testováno na Matlab 2020b.

% 1. Skript nacte anatomicky redukovane kovariaty ve forme 1D vektor zvoleneho subjekta.
% 2. Vektory jsou poskladany do spolecne matice.
% 3. Je vyuzita funkce corcoeff pro zskani matice Pearsonovych korelacnich koeficientu

%% Initialization and data loading, Inicializace a nacteni dat 
clc; clear all; close all; 

path='pac7.mat'; % ko1, pac2 ... (select subject, vyber subjekta)

front_data=load(['Anat_pca\Frontal\Pac\' path]);front_data=front_data.pwr;
pariet_data=load(['Anat_pca\Parietal\Pac\' path]);pariet_data=pariet_data.pwr;
temp_data=load(['Anat_pca\Temporal\Pac\' path]);temp_data=temp_data.pwr;
lefthem_data=load(['Anat_pca\Left_hem\Pac\' path]);lefthem_data=lefthem_data.pwr;
righthem_data=load(['Anat_pca\Right_hem\Pac\' path]);righthem_data=righthem_data.pwr;
occip_data=load(['Anat_pca\Occipital\Pac\' path]);occip_data=occip_data.pwr;

%% correlation analysis, korelacni analyza

pasma=[{'alfa'} {'beta'} {'delta'} {'theta'}];

% matrix for every subject, vytvoreni spolecne matice
front_mat=zeros(length(front_data.alfa),4);
temp_mat=zeros(length(temp_data.alfa),4);
occip_mat=zeros(length(occip_data.alfa),4);
lefthem_mat=zeros(length(lefthem_data.alfa),4);
pariet_mat=zeros(length(pariet_data.alfa),4);
righthem_mat=zeros(length(righthem_data.alfa),4);

% insertion, vlozeni do spolecne matice
for i=1:length(pasma)
front_mat(:,i)=front_data.(pasma{i});
temp_mat(:,i)=temp_data.(pasma{i});
occip_mat(:,i)=occip_data.(pasma{i});
lefthem_mat(:,i)=lefthem_data.(pasma{i});
righthem_mat(:,i)=righthem_data.(pasma{i});
pariet_mat(:,i)=pariet_data.(pasma{i});
end

% correlation analysis, Korelacni analyza 
mat=([front_mat temp_mat occip_mat lefthem_mat righthem_mat pariet_mat]); 
kor= corrcoef(mat); % kor je matice koeficientu pro dalsi analyzu


