%% EGI electrodes by anatomical regions%%
 
% Author:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Tested on Matlab 2020b.

% A set of commands for faster EGI 256 helmet channel selection.
%-------------------------------------------------------------------------%
% EGI elektrody dle anatomickych regionu

% Autor:  Bc. Filip Cerny, cernyfi4@fbmi.cvut.cz
% Testov�no na Matlab 2020b.

% soubor prikazu po rychle oznaceni anatomickych regionu z EGI256

el_map.right_hem={'E1','E2','E3','E4','E5','E6','E7','E10','E11','E12','E13',...
    'E14','E18','E19','E20','E25','E220','E221','E222','E223','E224','E215',...
    'E207','E198','E186','E185','E210','E211','E212','E213','E214','E206',...
    'E205','E204','E203','E202','E193','E192','E194','E195','E196','E197',...
    'E179','E180','E181','E182','E183','E184','E170','E171','E172','E173',...
    'E164','E160','E161','E162','E163','E155','E144','E132','E150','E151',...
    'E152','E153','E154','E143','E131','E138','E139','E140','E141','E142',...
    'E127','E128','E129','E130','E149','E159','E169','E178','E191','E31',...
    'E26','E21','E15','E8','E81','E90','E101','E119','E126','E137','E147'};

el_map.left_hem={'E32','E37','E46','E54','E61','E68','E83','E95','E106',...
    'E115','E124','E33','E38','E47','E55','E62','E69','E74','E84','E96',...
    'E107','E116','E125','E27','E34','E39','E48','E56','E63','E70','E75',...
    'E85','E97','E108','E117','E28','E35','E40','E49','E57','E64','E71',...
    'E76','E86','E98','E109','E118','E22','E29','E36','E41','E50','E58',...
    'E65','E72','E77','E87','E99','E110','E23','E30','E42','E51','E59',...
    'E66','E78','E88','E100','E16','E24','E43','E52','E60','E79','E89',...
    'E17','E44','E53','E80','E9','E45','E31','E26','E21','E15','E8','E81',...
    'E90','E101','E119','E126','E137','E147'};

el_map.occip={'E106','E107','E108','E113','E114','E115','E116','E117','E118',...
    'E120','E121','E122','E123','E124','E125','E126','E127','E133','E134',...
    'E135','E136','E137','E138','E139','E145','E146','E147','E148','E149',...
    'E150','E151','E156','E157','E158','E159','E160','E165','E166','E167',...
    'E168','E169','E174','E175','E176','E187'}; %S/N

el_map.left_temporal={'E111','E112','E105','E104','E103','E102','E92',...
    'E93','E94','E95','E83'};

el_map.right_temporal={'E208','E209','E199','E200','E201','E188','E189','E190',...
    'E191','E177','E178'};

el_map.temporal={'E111','E112','E105','E104','E103','E102','E92',...
    'E93','E94','E95','E83','E208','E209','E199','E200','E201','E188','E189','E190',...
    'E191','E177','E178'};

el_map.frontal={'E1','E2','E3','E4','E5','E6','E7','E8','E10','E11','E12',...
    'E13','E14','E15','E16','E17','E18','E19','E20','E21','E22','E23','E24',...
    'E25','E26','E27','E28','E29','E30','E31','E32','E33','E34','E35','E36',...
    'E37','E38','E39','E40','E41','E42','E43','E46','E47','E48','E49','E50',...
    'E51','E54','E55','E56','E57','E58','E61','E62','E63','E64','E68','E69'}; 

el_map.parietal={'E9','E44','E45','E52','E53','E59','E60','E65','E66','E70',...
    'E71','E72','E74','E75','E76','E77','E78','E79','E80','E81',...
    'E84','E85','E86','E87','E88','E89','E90','E96','E97','E98','E99','E100',...
    'E101','E109','E110','E119','E128','E129','E130','E131','E132','E140',...
    'E141','E142','E143','E144','E152','E153','E154','E155','E161','E162',...
    'E163','E164','E170','E171','E172','E173','E179','E180','E181','E182',...
    'E183','E184','E185','E186'};%  COM, REF

el_map.jaw={'E241','E242','E243','E244','E245','E246','E247','E248','E249',...
    'E250','E251','E252','E253','E254','E255','E256','E238','E239','E240',...
    'E234','E235','E236','E237','E230','E231','E232','E233','E255','E227','E228','E229'};


save elektrody_map.mat el_map;