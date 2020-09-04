%% SNF structural covariance analysis

addpath('../../data/BCT/2019_03_03_BCT')

location = "../../output/4_Out-of-model_features/Structural_Covariance/Measure_summaries/*.txt";

% import all residual files 
%% Getting global efficiency differences between each group

corr_1 = corrcoef(group1resid);
corr_2 = corrcoef(group2resid);
corr_3 = corrcoef(group3resid);
corr_4 = corrcoef(group4resid);

corr_adhd = corrcoef(ADHDresid);
corr_asd = corrcoef(ASDresid);
corr_ocd = corrcoef(OCDresid);

% Creating empty tables
group_diff=zeros(6,5);
group_eff=zeros(4,5);
clin_diff=zeros(6,3);
clin_eff=zeros(3,3);
group_stren_diff=zeros(6,5);
group_stren=zeros(4,5);
clin_stren_diff=zeros(6,3);
clin_stren=zeros(4,5);

 

% Calculating global efficiency and network strength for data-driven and dx groups at 9 
%% different thresholds between 10 and 50 percent of network edges

row=0
for threshold=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    threshold
    row=row+1
    
    % thresholding correlation matrices
    W_1 = threshold_proportional(corr_1, threshold);
    W_2 = threshold_proportional(corr_2, threshold);
    W_3 = threshold_proportional(corr_3, threshold);
    W_4 = threshold_proportional(corr_4, threshold);
    
    C_adhd = threshold_proportional(corr_adhd, threshold);
    C_asd = threshold_proportional(corr_asd, threshold);
    C_ocd = threshold_proportional(corr_ocd, threshold);
    
    %calculating the mean network strength for each group
    
    stren_1=triu(W_1, 1); %getting the upper part of the triangle
    stren_1=stren_1(:); %organizing into a column
    stren_1(any(~stren_1,2),:) = []; %removing all that aren't zeros
    S_1=mean(stren_1); %calculating the average
    
    stren_2=triu(W_2, 1);
    stren_2=stren_2(:);
    stren_2(any(~stren_2,2),:) = [];
    S_2=mean(stren_2);
    
    stren_3=triu(W_3, 1);
    stren_3=stren_3(:);
    stren_3(any(~stren_3,2),:) = [];
    S_3=mean(stren_3);
    
    stren_4=triu(W_4, 1);
    stren_4=stren_4(:);
    stren_4(any(~stren_4,2),:) = [];
    S_4=mean(stren_4);
    
    stren_adhd=triu(C_adhd, 1);
    stren_adhd=stren_adhd(:);
    stren_adhd(any(~stren_adhd,2),:) = [];
    S_adhd=mean(stren_adhd);
    
    stren_asd=triu(C_asd, 1);
    stren_asd=stren_asd(:);
    stren_asd(any(~stren_asd,2),:) = [];
    S_asd=mean(stren_asd);
    
    stren_ocd=triu(C_ocd, 1);
    stren_ocd=stren_ocd(:);
    stren_ocd(any(~stren_ocd,2),:) = [];
    S_ocd=mean(stren_ocd);
       
    %binarizing the networks once the average network strength has been calculated
    W_1(W_1>0)=1;
    W_2(W_2>0)=1;
    W_3(W_3>0)=1;
    W_4(W_4>0)=1;
    
    C_adhd(C_adhd>0)=1;
    C_asd(C_asd>0)=1;
    C_ocd(C_ocd>0)=1;

    % calculating global efficiency 
    E_1=efficiency_bin(W_1,0) %0 for global, 1 for local
    E_2=efficiency_bin(W_2,0)
    E_3=efficiency_bin(W_3,0)
    E_4=efficiency_bin(W_4,0)
    
    EC_adhd=efficiency_bin(C_adhd,0) %0 for global, 1 for local
    EC_asd=efficiency_bin(C_asd,0)
    EC_ocd=efficiency_bin(C_ocd,0)
    
    % rounding the efficiency measure to 3 digits
    E_1=round(E_1, 3); 
    E_2=round(E_2, 3);
    E_3=round(E_3, 3);
    E_4=round(E_4, 3);
    
    EC_adhd=round(EC_adhd, 3); 
    EC_asd=round(EC_asd, 3);
    EC_ocd=round(EC_ocd, 3);

    % calculating differences in efficiency between groups
    group_diff(row,1)=abs(E_1 - E_2);
    group_diff(row,2)=abs(E_1 - E_3);
    group_diff(row,3)=abs(E_1 - E_4);
    group_diff(row,4)=abs(E_3 - E_2);
    group_diff(row,5)=abs(E_4 - E_2);
    group_diff(row,6)=abs(E_3 - E_4);
    
    clin_diff(row,1)=abs(EC_adhd - EC_asd);
    clin_diff(row,2)=abs(EC_adhd - EC_ocd);
    clin_diff(row,3)=abs(EC_asd - EC_ocd);
    
    % calculating mean network strength differences between the groups
    group_stren_diff(row,1)=abs(S_1 - S_2);
    group_stren_diff(row,2)=abs(S_1 - S_3);
    group_stren_diff(row,3)=abs(S_1 - S_4);
    group_stren_diff(row,4)=abs(S_3 - S_2);
    group_stren_diff(row,5)=abs(S_4 - S_2);
    group_stren_diff(row,6)=abs(S_3 - S_4);
    
    clin_stren_diff(row,1)=abs(S_adhd - S_asd);
    clin_stren_diff(row,2)=abs(S_adhd - S_ocd);
    clin_stren_diff(row,3)=abs(S_ocd - S_asd);
    
    % filling in the tables for this threshold
    group_eff(row, 1)=E_1;
    group_eff(row, 2)=E_2;
    group_eff(row, 3)=E_3;
    group_eff(row, 4)=E_4;
    
    group_stren(row, 1)=S_1;
    group_stren(row, 2)=S_2;
    group_stren(row, 3)=S_3;
    group_stren(row, 4)=S_4;
    
    clin_eff(row, 1)=EC_adhd;
    clin_eff(row, 2)=EC_asd;
    clin_eff(row, 3)=EC_ocd;
    
    clin_stren(row, 1)=S_adhd;
    clin_stren(row, 2)=S_asd;
    clin_stren(row, 3)=S_ocd;
end

% Creating tables summarizing global efficiencies for groups
col_names={'diff_12', 'diff_13', 'diff_14', 'diff_23', 'diff_24', 'diff_34'}
row_names={'0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.45', '0.5'}
group_GE_diffs= array2table(group_diff,'RowNames',row_names,'VariableNames',col_names)

col_names={'group1', 'group2', 'group3', 'group4', 'blank'}
group_GEs= array2table(group_eff,'RowNames',row_names,'VariableNames',col_names)

% Creating tables summarizing network strength for groups
col_names={'diff_12', 'diff_13', 'diff_14', 'diff_23', 'diff_24', 'diff_34'}
group_stren_diff= array2table(group_stren_diff,'RowNames',row_names,'VariableNames',col_names)

col_names={'group1', 'group2', 'group3', 'group4', 'blank'}
group_strengths= array2table(group_stren,'RowNames',row_names,'VariableNames',col_names)

% Diagnostic groups - global efficiency
col_names={'diff_adhd_asd', 'diff_adhd_ocd', 'diff_asd_ocd'}
clin_GE_diffs= array2table(clin_diff,'RowNames',row_names,'VariableNames',col_names)

col_names={'adhd', 'asd', 'ocd'}
clin_GE= array2table(clin_eff,'RowNames',row_names,'VariableNames',col_names)

% Diagnostic groups - network strength
col_names={'diff_adhd_asd', 'diff_adhd_ocd', 'diff_asd_ocd'}
clin_strength_diff= array2table(clin_stren_diff,'RowNames',row_names,'VariableNames',col_names)

col_names={'adhd', 'asd', 'ocd', 'blank1', 'blank2'}
clin_strengths= array2table(clin_stren, 'RowNames',row_names,'VariableNames',col_names)


%% writing out global efficiency and network strength tables
writetable(location, group_GE_diffs)
writetable(location, group_GEs)

writetable(location, group_stren_diff) 
writetable(location, group_strengths)

writetable(location, clin_GE_diffs)
writetable(location, clin_GE)

writetable(location, clin_strength_diff)
writetable(location, clin_strengths)

%% Getting network densities and differences between groups

% setting up empty tables
group_density_diff=zeros(6,6);
group_density=zeros(4,4);
clin_density_diff=zeros(6,3);
clin_density=zeros(3,3);

% thresholding networks at a range of absolute network strength thresholds
row=0
for threshold=0.33:0.01:0.55
    threshold
    row=row+1

% setting up structural covariance networks
corr_1 = corrcoef(group1resid);
corr_2 = corrcoef(group2resid);
corr_3 = corrcoef(group3resid);
corr_4 = corrcoef(group4resid);

corr_adhd = corrcoef(ADHDresid);
corr_asd = corrcoef(ASDresid);
corr_ocd = corrcoef(OCDresid);

    % absolute thresholding correlation matrices
    % overwrite W_1 each time
    W_1 = threshold_absolute(corr_adhd, threshold);
    W_1(W_1>0)=1; % binarizing them
    D_adhd=density_und(W_1); %calculating density
    
    W_1 = threshold_absolute(corr_asd, threshold);
    W_1(W_1>0)=1; % binarizing them
    D_asd=density_und(W_1);
    
    W_1 = threshold_absolute(corr_ocd, threshold);
    W_1(W_1>0)=1; % binarizing them
    D_ocd=density_und(W_1);

    W_1 = threshold_absolute(corr_1, threshold);
    W_1(W_1>0)=1; % binarizing them
    D_1=density_und(W_1);
    
    W_2 = threshold_absolute(corr_2, threshold);
    W_2(W_2>0)=1; % binarizing them
    D_2=density_und(W_2);
    
    W_3 = threshold_absolute(corr_3, threshold);
    W_3(W_3>0)=1; % binarizing them
    D_3=density_und(W_3);
    
    W_4 = threshold_absolute(corr_4, threshold);
    W_4(W_4>0)=1; % binarizing them
    D_4=density_und(W_4);
    
    % filling in tables with values
    clin_density(row, 1)=D_adhd;
    clin_density(row, 2)=D_asd;
    clin_density(row, 3)=D_ocd;
    
    % calculating density differences between groups
    clin_density_diff(row,1)=abs(D_adhd - D_asd);
    clin_density_diff(row,2)=abs(D_adhd - D_ocd);
    clin_density_diff(row,3)=abs(D_asd - D_ocd);

    group_density(row, 1)=D_1;
    group_density(row, 2)=D_2;
    group_density(row, 3)=D_3;
    group_density(row, 4)=D_4;
    
    group_density_diff(row,1)=abs(D_1 - D_2);
    group_density_diff(row,2)=abs(D_1 - D_3);
    group_density_diff(row,3)=abs(D_1 - D_4);
    group_density_diff(row,4)=abs(D_3 - D_2);
    group_density_diff(row,5)=abs(D_4 - D_2);
    group_density_diff(row,6)=abs(D_3 - D_4);
    
end

col_names={'G1', 'G2', 'G3', 'G4'}
group_density=array2table(group_density, 'VariableNames',col_names)

col_names={'diff_12', 'diff_13', 'diff_14', 'diff_32', 'diff_42', 'diff_34'}
group_density_diffs= array2table(group_density_diff, 'VariableNames',col_names)

col_names={'adhd', 'asd', 'ocd'}
clin_density= array2table(clin_density, 'VariableNames',col_names)

col_names={'adhd_asd', 'adhd_ocd', 'asd_ocd'}
clin_density_diffs= array2table(clin_density_diff, 'VariableNames',col_names)

writetable(location, group_density)
writetable(location, group_density_diffs)
writetable(location, clin_density)
writetable(location, clin_density_diffs)

%% Permutation testing to determine if differences in global efficiency and mean network
% strength are significant between data-driven groups

% setting up empty tables
eff_table=zeros(4,4);
ge = zeros(1000, 1);
stren_table=zeros(4,4);
stren = zeros(1000, 1);

% randomizing each cortical thickness matrix 1000 times, then creating a structural 
% covariance matrix, then calculating global efficiency and
% finding group differences to create a null distribution

all=[group1resid; group2resid; group3resid; group4resid];


col=0


for group_1={17, 47, 31, 41} % data-driven group sizes
    g1_size=group_1{1}
    
    
for group_2={17, 47, 31, 41}
    g2_size=group_2{1}
    col=col+1
    row=0;
    
for threshold=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5] % percent of edges
    threshold
    row=row+1
 
for i=1:1000
    idx=randperm(136); %randomly ordering 1-136 participants
    rand_mat=all(idx,:);%rearranging to match random ordering    
       
    g1=rand_mat(1:g1_size,:);
    g2=rand_mat(g1_size+1:g1_size+g2_size,:);


    corr_1 = corrcoef(g1); %creating structural covariance matrix
    W_1 = threshold_proportional(corr_1, threshold); %thresholding
    W_1(W_1>0)=1; %binarizing
    E_1=efficiency_bin(W_1,0); %calculating efficiency

    corr_1=triu(corr_1, 1); %taking the top of the triangle
    corr_1=corr_1(:); %turning into a single column
    corr_1(any(~corr_1,2),:) = []; %removing zeros
    S_1=mean(corr_1); %calculating mean strength
    
    corr_2 = corrcoef(g2);
    W_2 = threshold_proportional(corr_2, threshold);
    W_2(W_2>0)=1;
    E_2=efficiency_bin(W_2,0);
    corr_2=triu(corr_2, 1);
    corr_2=corr_2(:);
    corr_2(any(~corr_2,2),:) = [];
    S_2=mean(corr_2);
    
    diff = E_1 - E_2; % calculating difference in efficiency
    diff = abs(diff); % absolute valuing difference
    
    S_diff = abs(S_1 - S_2); %calculating difference in strenth
    
    % recording differences
    ge(i,1)=diff; 
    stren(i, 1)=S_diff;
    
end
    % across permutations, calculating value of significant difference p=0.05
    threshold_100 = quantile(ge, 0.95);
    threshold_stren = quantile(stren, 0.95);
    % recording significant differences
    eff_table(row,col)=threshold_100
    stren_table(row,col)=threshold_stren
end

end

end

col_names={'group1_1', 'group2_1', 'group3_1', 'group4_1', 'group1_2', 'group2_2', 'group3_2', 'group4_2', 'group1_3', 'group2_3', 'group3_3', 'group4_3', 'group1_4', 'group2_4', 'group3_4', 'group4_4'}
row_names={'0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.45', '0.5'}
group_sig_ge_diff_1000_95 = array2table(eff_table,'RowNames',row_names,'VariableNames',col_names)
group_sig_stren_diff_1000_95 = array2table(stren_table,'RowNames',row_names,'VariableNames',col_names)

writetable(location, group_sig_stren_diff_1000_95)

writetable(location, group_sig_ge_diff_1000_95)


%% Permutation testing to determine significant density differences between 
%data-driven groups

% setting up empty tables
density_table=zeros(4,4);
ge = zeros(1000, 1);

all=[group1resid; group2resid; group3resid; group4resid]


col=0
for group_1={17, 47, 31, 41} % data-driven group sizes
    g1_size=group_1{1}
    
    
for group_2={17, 47, 31, 41}
    g2_size=group_2{1}
    col=col+1
    row=0;
    
for threshold=0.34:0.01:0.55
    threshold
    row=row+1
 
for i=1:1000
    idx=randperm(136); %randomly reordering participants
    rand_mat=all(idx,:);    
       
    g1=rand_mat(1:g1_size,:); %creating two cortical thickness matrices 
    g2=rand_mat(g1_size+1:g1_size+g2_size,:);

    g_cor_1 = corrcoef(g1); %creating structural covariance matrix
    ab_1 = threshold_absolute(g_cor_1, threshold);
    ab_1(ab_1>0)=1; %binarizing matrices
    den_1=density_und(ab_1); %calculating network density
        
    g_cor_2 = corrcoef(g2);
    ab_2 = threshold_absolute(g_cor_2, threshold);
    ab_2(ab_2>0)=1;
    den_2=density_und(ab_2);
        
    diff = den_1 - den_2; %calculating difference between groups
    diff = abs(diff);

    ge(i,1)=diff; %recording difference
    
end
    threshold_density = quantile(ge, 0.95); %calculating significant difference
    density_table(row,col)=threshold_density %recording significant difference
end

end

end


col_names={'group1_1', 'group2_1', 'group3_1', 'group4_1', 'group1_2', 'group2_2', 'group3_2', 'group4_2', 'group1_3', 'group2_3', 'group3_3', 'group4_3', 'group1_4', 'group2_4', 'group3_4', 'group4_4', 'blank'}
group_sig_density_1000_95 = array2table(density_table)

writetable(location, group_sig_density_1000_95)


%% Permutation testing to determine significant density differences between 
% diagnostic groups

eff_table=zeros(9,6); %creating empty tables
ge = zeros(1000, 1);
stren_table=zeros(9,6);
stren = zeros(1000, 1);

% randomizing each cortical thickness matrix 1000 times, then creating a structural 
% covariance matrix, then calculating global efficiency and
% finding group differences to create a null distribution

all=[ADHDresid; ASDresid; OCDresid];


col=0

for group_1={47, 72, 17} %diagnostic group sizes
    g1_size=group_1{1}
 
for group_2={47, 17}
    g2_size=group_2{1}
    col=col+1
    row=0
    
for threshold=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    threshold
    row=row+1
 
for i=1:1000
    idx=randperm(136); % randomly shuffling numbers 
    rand_mat=all(idx,:); % reordering the matrix of participants
    
    g1=rand_mat(1:g1_size,:);
    g2=rand_mat(g1_size+1:g1_size+g2_size,:);

    corr_s_1 = corr(g1); %creating structural covariance matrix
    W_1 = threshold_proportional(corr_s_1, threshold); %thresholding
    W_1(W_1>0)=1; %binarizing
    E_1=efficiency_bin(W_1,0); %calculating global efficiency
    corr_s_1=triu(corr_s_1, 1); %getting the top triangle of the matrix
    corr_s_1=corr_s_1(:); %creating one column
    corr_s_1(any(~corr_s_1,2),:) = []; %removing zeros
    S_1=mean(corr_s_1); % calculating mean network strength
    
    corr_s_2 = corr(g2);
    W_2 = threshold_proportional(corr_s_2, threshold);
    W_2(W_2>0)=1;
    E_2=efficiency_bin(W_2,0); % global efficiency
    corr_s_2=triu(corr_s_2, 1);
    corr_s_2=corr_s_2(:);
    corr_s_2(any(~corr_s_2,2),:) = [];
    S_2=mean(corr_s_2); %mean strength
    %imagesc(W_1)
    
    diff = E_1 - E_2; %calculating difference between groups in efficiency
    diff = abs(diff); %absolute value

    S_diff = abs(S_1 - S_2); % calculating difference in network strength

    ge(i,1)=diff; %recording differences
    stren(i, 1)=S_diff;
    
end
    threshold_ge = quantile(ge, 0.95); %calculating significance threshold
    threshold_stren = quantile(stren, 0.95);
    eff_table(row,col)=threshold_ge % recording significant threshold for edge threshold
    stren_table(row,col)=threshold_stren
end

end

end

col_names={'adhd_adhd', 'adhd_ocd', 'asd_adhd', 'asd_ocd', 'ocd_adhd', 'ocd_ocd'}
row_names={'0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.45', '0.5'}
clinical_ge_diff_1000_95 = array2table(eff_table,'RowNames',row_names,'VariableNames',col_names)
clinical_strength_diff_1000_95 = array2table(stren_table,'RowNames',row_names,'VariableNames',col_names)


writetable(location, clinical_ge_diff_1000_95)
writetable(location, clinical_strength_diff_1000_95)

%% Calculating significant density difference between diagnostic groups

density_table=zeros(4,6);
ge = zeros(1000, 1);

% randomizing each cortical thickness matrix 1000 times, then creating a structural 
% covariance matrix, then calculating density and
% finding group differences to create a null distribution

all=[ADHDresid; ASDresid; OCDresid];


col=0
for group_1={41, 78, 17}
    g1_size=group_1{1}
 
for group_2={41, 17}
    g2_size=group_2{1}
    col=col+1
    row=0
    
for threshold=0.33:0.01:0.55
    threshold
    row=row+1
 
for i=1:1000
    idx=randperm(136);
    rand_mat=all(idx,:);    
    %idx=randperm(68);
    %rand_mat=rand_mat(:,idx);
       
    g1=rand_mat(1:g1_size,:);
    g2=rand_mat(g1_size+1:g1_size+g2_size,:);

    g_cor_1 = corr(g1); %creating matrix
    %corr_1 = atanh(corr_1);
    ab_1 = threshold_absolute(g_cor_1, threshold);
    ab_1(ab_1>0)=1;
    den_1=density_und(ab_1);
    
    %imagesc(W_1)
    
    g_cor_2 = corr(g2);
    ab_2 = threshold_absolute(g_cor_2, threshold);
    ab_2(ab_2>0)=1;
    den_2=density_und(ab_2);
        
    diff = den_1 - den_2;
    diff = abs(diff);

    ge(i,1)=diff;
    
    
end
    threshold_100 = quantile(ge, 0.95);
    density_table(row,col)=threshold_100
end

end

end

col_names={'adhd_adhd', 'adhd_ocd', 'asd_adhd', 'asd_ocd', 'ocd_adhd', 'ocd_ocd'}
clinical_density_diff_1000_95 = array2table(density_table,'VariableNames',col_names)

writetable(location, clinical_density_diff_1000_95)





