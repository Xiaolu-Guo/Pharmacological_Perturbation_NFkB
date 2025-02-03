clear
stim = {'TNF','LPS','CpG','polyIC','Pam3CSK'};
dose = [0.3 3.3 33; 0.33 3.3 33; 33 100 330; 3.3 33 100; 10 33 100];
dose_scale = [1/5200, 1/24000, 1/1000, 1000/5e6, 1/1500];

names{1} = {'TNF','TNFR','TNFR_TNF','C1','IkBat','IkBa','IkBan','IKK','NFkBn','TAK1', 'IkBaNFkBn', 'NFkB','IkBaNFkB'};
names{2} = {'TLR4','TLR4LPS','TLR4LPSen','TRIF','MyD88','TRAF6','IKK','IkBat','NFkBn','TAK1','IkBaNFkBn', 'NFkB','IkBaNFkB'};
names{3} = {'CpG','CpG_en','TLR9','TLR9_CpG','TLR9_N','MyD88','TRAF6','IKK','NFkBn','TAK1','IkBaNFkBn', 'NFkB','IkBaNFkB'};
names{4} = {'polyIC','polyIC_en','TLR3','TLR3_polyIC','TRIF','TRAF6','IKK','NFkBn','TAK1','IkBaNFkBn', 'NFkB','IkBaNFkB'};
names{5} = {'CD14',  'Pam3CSK', 'CD14_P3CSK','TLR2','TLR2_P3CSK','MyD88','TRAF6','TAK1','IKK','NFkBn','IkBaNFkBn', 'NFkB','IkBaNFkB'};
%species2save = {'NFkBn','IkBaNFkBn', 'NFkB', 'IkBaNFkB', 'IkBa', 'IkBan'};
species2save = {'NFkBn','IkBaNFkBn'}; 

sa_choice = 2; %1 for biochem
sa_type = {'Biochem', 'Pharm'};

options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

opts1 = detectImportOptions('NFkB_param_groups_EB.xlsx','Sheet','Param_Groupings');

for jj = 1: length(opts1.VariableNames)
    opts1 = setvartype(opts1, opts1.VariableNames{jj}, 'char');
end

params_table = readtable('NFkB_param_groups_EB.xlsx',opts1);

rxn_num = str2double(params_table.Reaction_Number(1:94));
param_num = params_table.Param_Number(1:94);
if sa_choice == 1
    biochem_group = params_table.Biochem_Group(1:94);
    [unique_bc_group, ~, ubcg_ind] = unique(biochem_group);
    mparams = 10.^linspace(-1, 1, 21); %multiply parameters by 0.1-10X
elseif sa_choice == 2
    biochem_group = params_table.DID(1:94);
    [unique_bc_group, ~, ubcg_ind] = unique(biochem_group);
    unique_bc_group(strcmp('',unique_bc_group)) = [];
    %ubcg_ind(find(ubcg_ind == 1)) = [];
    mparams = 10.^linspace(-3, 0, 21); %multiply parameters by 0.001-1X
end

%NOTE: if doing pharm sens analysis, indexing starts at 2!

bc_sens_analysis = struct;
c = 0;
for bcg_ind = 1:length(unique_bc_group)
    
    if sa_choice == 1
        rxn2mod_ind = find(ubcg_ind == bcg_ind);
        rxn2mod = unique(rxn_num(rxn2mod_ind), 'stable');
    elseif sa_choice == 2 && bcg_ind <= 10
        rxn2mod_ind = find(ubcg_ind == (bcg_ind + 1));
        rxn2mod = unique(rxn_num(rxn2mod_ind), 'stable');
        if bcg_ind == 10 %TSA
            mparams = 10.^linspace(0, 3, 21);
        else
            mparams = 10.^linspace(-3, 0, 21);
        end
    end
    
    for mp_ind = 1:length(mparams)
        
        [v.PARAMS, v.SPECIES] = nfkbInitialize();
        options.v.SPECIES = v.SPECIES;
        %options0.v.SPECIES = v.SPECIES;
        %options0.v.PARAMS = v.PARAMS;
        v.PARAMS(rxn2mod,1) = v.PARAMS(rxn2mod,1) * mparams(mp_ind);
        options.v.PARAMS = v.PARAMS;
        
        for stim_ind = 1:length(stim)
            output = [];
            for dose_ind = 1:length(dose(1,:))
                
                % Simulate all doses (only need to equilibrate on first iteration)
                if isempty(output)
                   [t,x,simdata] = nfkbSimulate({stim(stim_ind),dose(stim_ind,dose_ind)*dose_scale(stim_ind)},...
                                   names{stim_ind},[], {},options);
                else
                    options.STEADY_STATE = simdata.STEADY_STATE;
                    [~,x] = nfkbSimulate({stim(stim_ind),dose(stim_ind,dose_ind)*dose_scale(stim_ind)},...
                            names{stim_ind},[], {},options);
                end
                output = cat(3,output,x);
                disp([unique_bc_group(bcg_ind)+"X"+num2str(mparams(mp_ind))+...
                    "; "+stim(stim_ind)+" at dose "+ dose(stim_ind,dose_ind)*dose_scale(stim_ind)]);
                %clear options.STEADY_STATE
            end
            for species_ind = 1:length(species2save)
                for dose_ind2 = 1:length(dose(1,:))
                    c = c + 1;
                    bc_sens_analysis.bc_group{c} = unique_bc_group{bcg_ind};
                    bc_sens_analysis.rxn_num{c} = rxn2mod;
                    bc_sens_analysis.param_num{c} = 1;
                    bc_sens_analysis.multiplier{c} = mparams(mp_ind);
                    bc_sens_analysis.multiplier_count{c} = mp_ind;
                    bc_sens_analysis.param_value{c} = v.PARAMS(rxn2mod,1);
                    bc_sens_analysis.ligand{c} = stim{stim_ind};
                    bc_sens_analysis.dose{c} = dose(stim_ind,dose_ind2);
                    bc_sens_analysis.dose_str{c} = num2str(dose(stim_ind,dose_ind2));
                    bc_sens_analysis.species{c} = species2save{species_ind};
                    bc_sens_analysis.trajectory{c} =  output(:,strcmp(names{stim_ind},species2save(species_ind)),dose_ind2)';
                end
                
            end
        end
        %Remove comment if want params to be changed for steady state calcs
        %TLR3 (PolyIC) doesn't converge at high doses of chlorquine...
        %now introduce drug before steady state
        options = rmfield(options, 'STEADY_STATE');
    end
end


bc_sens_table = table(bc_sens_analysis.bc_group(:),bc_sens_analysis.rxn_num(:),bc_sens_analysis.param_num(:),...
                bc_sens_analysis.multiplier(:),bc_sens_analysis.multiplier_count(:),bc_sens_analysis.param_value(:),bc_sens_analysis.ligand(:),...
                bc_sens_analysis.dose(:),bc_sens_analysis.dose_str(:),bc_sens_analysis.species(:),bc_sens_analysis.trajectory(:),'VariableNames',fieldnames(bc_sens_analysis));

            
         
data_ttl = strcat(sa_type{sa_choice},'_all_stim_0.001-1X_simulations.mat');
save(data_ttl,'bc_sens_analysis','bc_sens_table');   