function biodomain = biochemical_domains(protein)

% Code defining the experimentally known biochemical domains in each
% protein
% 
% Last updated: 2017-03-25
% 
% INPUT
%     protein               Options:'Gag','Nef','NS34A','NS4B'

%%
if strcmp(protein,'Gag') 
    
    %% HIV Gag

    
    load('sites_involved_in_interface_hexamer_p24_7.mat') %PDB:3GV2
    gag_p24_intra_hex_int = sites_involved_in_interface_hexamer+132;
    load('sites_involved_in_interface_pentamer_p24_7.mat') %PDB:3P05
    gag_p24_intra_pent_int = sites_involved_in_interface_pentamer+132;
    load('sites_involved_in_interface_between_hexamers_p24_7.mat') %PDB:2KOD
    gag_p24_inter_hex_int = sites_involved_in_interface_between_hexamers+132;
    
    %Epitopes targeted by elite controllers
    gag_epitopes_elite = unique([145:155 203:212 240:249 162:172 147:155 263:272 298:306 19:27 180:188 305:313]); %Supplementary Table 3
    %Epitopes targeted by progressors
    gag_epitopes_progressors = unique([78:86 355:363 260:267 254:262 36:44]); %Supplementary Table 3
    
    biodomain(1).name = 'P17-Mem-Bin-Dom';
    biodomain(1).sites = 1:31;
    biodomain(2).name = 'P24-SP1-Int';
    biodomain(2).sites = 360:368;
    biodomain(3).name = 'P24-Intra-Hex-Int';
    biodomain(3).sites = gag_p24_intra_hex_int;
    biodomain(4).name = 'P7-Zinc-Fingers';
    biodomain(4).sites = [392 395 400 405 406:412 413 416 421 426];
    biodomain(5).name = 'P24-Inter-Hex-Int';
    biodomain(5).sites = gag_p24_inter_hex_int;
    biodomain(6).name = 'Epitopes associated with viral contol';
    biodomain(6).sites = gag_epitopes_elite;
    biodomain(7).name = 'Epitopes associated with disease progression';
    biodomain(7).sites = gag_epitopes_progressors;
    
elseif strcmp(protein,'Nef')
    
    %% HIV Nef
    
    biodomain(1).name = 'Enh-Vir-Inf';
    biodomain(1).sites = 69:78; 
    biodomain(2).name = 'HLA1-Down-Reg';
    biodomain(2).sites = [62:65 72 75 77 82 86]; 
    biodomain(3).name = 'Intra-Dim-Int';
    biodomain(3).sites = [105, 108, 109, 112, 115, 116, 121, 122, 123]; 
    biodomain(4).name = 'CD4-Down-Reg';
    biodomain(4).sites = [57:59 95:97 106 109 110]; 
    
elseif strcmp(protein,'NS34A')
    
    %% HCV NS3-4A
    
    id_chainA = [435:453 477:488 524:536];
    id_chainB = [545:553 584:591];
    
    biodomain(1).name = 'NS3-NS4A-Pro-Act';
    biodomain(1).sites = [1:22 652:665]; 
    biodomain(2).name = 'NS3-NS4A-Mem-Asso';
    biodomain(2).sites = [10:24 [1:20]+631]; 
    biodomain(3).name = 'NS5A-Hyper-Phos';
    biodomain(3).sites = unique([1027:1036 1165 1678 1690:1692 1697:1699 1703:1705 1707:1710]-1026); 
    biodomain(4).name = 'NS3-Motif-Enz-Heli';
    biodomain(4).sites = [1486:1493]-1026;
    biodomain(5).name = 'NS3-Intra-Dimer-Int';
    biodomain(5).sites = [id_chainA id_chainB];
    
elseif strcmp(protein,'NS4B')
    
    %% HCV NS4B
    
    ns4b_H1 = 201:212; 
    ns4b_H2 = 229:253; 
    ns4b_bzip = 20:55;  
    ns4b_oligomer = [1:70 257 261];
    
    biodomain(1).name = 'Viral-Rep-Assm';
    biodomain(1).sites = ns4b_H1; 
    biodomain(2).name = 'NS4B-ATF6beta';
    biodomain(2).sites = ns4b_bzip; 
    biodomain(3).name = 'Oligomerization';
    biodomain(3).sites = ns4b_oligomer;
    
end