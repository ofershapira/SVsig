% load or generate ICGC model and data structures
model_exist = false;

WorkDir = pwd;
addpath(genpath(pwd))
%DataDir = strcat(WorkDir,'/data/');
%TracksDir = strcat(WorkDir,'/tracks/');

% load data table with merged SV with the following columns:
% {seqnames, start, strand1, altchr, altpos, strand2, subtype(histology), sv_id, sid(sample ID), donor_unique_id} 
sv_file=strcat(WorkDir,'/tracks/merged_1.6.1.csv');
SVTable=readtable(sv_file, 'Delimiter', '\t');

if model_exist

    load ICGC_2D_SV_model.mat

else

    mix_model

end


EventLengthThreshold=1e2;
len_filter=1e6;
bks_cluster=1;

% events array
events0=zeros(height(SVTable),12);

if sum(strcmp(SVTable.Properties.VariableNames, 'histology_abbreviation'))>0
    SVTable.Properties.VariableNames{'histology_abbreviation'} = 'subtype';
end

if isa(SVTable.seqnames,'numeric')    
    events0(:,1)=SVTable.seqnames;
else
    events0(:,1)=chr_str2num(SVTable.seqnames)';
end

if isa(SVTable.altchr,'numeric')
    events0(:,4)=SVTable.altchr;    
else
    events0(:,4)=chr_str2num(SVTable.altchr)';
end

events0(:,2)=SVTable.start;
[Ustrand1, ia_strand1, ic_strand1]=unique(SVTable.strand);
events0(:,3)=ic_strand1;


events0(:,5)=SVTable.altpos;
[Ustrand2, ia_strand2, ic_strand2]=unique(SVTable.altstrand);
events0(:,6)=ic_strand2;

[UTumor, ia_code, ic_code]=unique(SVTable.subtype);
events0(:,7)=ic_code;

[Uevent, ia_event, ic_event]=unique(SVTable.sv_id);
events0(:,8)=ic_event;

[Usample, ia_sample, ic_sample]=unique(SVTable.sid);
events0(:,9)=ic_sample;

if sum(strcmp('donor_unique_id',SVTable.Properties.VariableNames))>0
    [Upatient, ia_patient, ic_patient]=unique(SVTable.donor_unique_id);
else
    [Upatient, ia_patient, ic_patient]=unique(SVTable.sid);
end
events0(:,10)=ic_patient;

%events0(:,11)=SVTable.HOMLEN;
%events0(:,12)=SVTable.INSLEN;

disp(strcat('total events from vcfs: ',num2str(length(events0))));
% filter mask track
[events0,masked_events] = mask_events( events0,mask_track );
disp(strcat('total events after masked regions: ',num2str(length(events0))));

% events matrix
mfull0=sparse(length(bins),length(bins));
bins_event_tble0=zeros(length(events0),3);
for c1 = 1:length(events0),
    bini0=find(bins(:,1)==events0(c1,1) & bins(:,2)<=events0(c1,2) & bins(:,3)>=events0(c1,2));
    binj0=find(bins(:,1)==events0(c1,4) & bins(:,2)<=events0(c1,5) & bins(:,3)>=events0(c1,5));
    if ~isempty(bini0) && ~isempty(binj0)
        bins_event_tble0(c1,1) = bini0;
        bins_event_tble0(c1,2) = binj0;
        bini = min(bins_event_tble0(c1,1),bins_event_tble0(c1,2));
        binj = max(bins_event_tble0(c1,1),bins_event_tble0(c1,2));
        bins_event_tble0(c1,3) = sub2ind([length(bins) length(bins)],bini,binj); % note: assign values only to upper tria of the matrix
        mfull0(bini,binj) = mfull0(bini,binj) + 1;
    end
end

% remove unassigned events
unassigned_events=bins_event_tble0(:,1)==0;
events0(unassigned_events,:)=[];
bins_event_tble0(unassigned_events,:)=[];
disp(strcat('total events after removing unassigned events: ',num2str(length(events0))));

% remove same bin same sample events
T_bin_sample=table();
T_bin_sample.bin=bins_event_tble0(:,3);
T_bin_sample.sample=events0(:,9);
[u_t,ia_t,ic_t]=unique(T_bin_sample);
events=events0(ia_t,:);
bins_event_tble=bins_event_tble0(ia_t,:);
disp(strcat('total events after removing unassigned events: ',num2str(length(events))));

% remove below length threshold events
below_length_th=(abs(events(:,5)-events(:,2))<EventLengthThreshold)&(events(:,1)==events(:,4));
events(below_length_th,:)=[];
bins_event_tble(below_length_th,:)=[];
disp(strcat('total events after removing below length threshold: ',num2str(length(events))));



% update mfull
mfull=sparse(length(bins),length(bins));
for c1 = 1:length(events),
    bini0=find(bins(:,1)==events(c1,1) & bins(:,2)<=events(c1,2) & bins(:,3)>=events(c1,2));
    binj0=find(bins(:,1)==events(c1,4) & bins(:,2)<=events(c1,5) & bins(:,3)>=events(c1,5));
    if ~isempty(bini0) && ~isempty(binj0)
        mfull(bini0,binj0) = mfull(bini0,binj0) + 1;
    end
end

if ~bks_cluster
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull+mfull', mix_model, [], [],1);
else
    sij1dx = length_dist_1d_bins(events,chsize,10);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull+mfull', mix_model, bins, events, sij1dx, chsize, CHR, 1);
end

[hitstable_mix,hitstable_mix_lookup] = HitsTableCV(mfull_pval_mix,pa_mix, pval_tophits_mix, bins_event_tble, qFDR_mix, events, refgene_tble);

CuratedFusionGene0=CuratedFusionGene(1:end-3,:);
TbyGene_mix = TophitsByGenes(hitstable_mix,hitstable_mix_lookup,1e4,bins,refgene,refgene_tble,[],CosmicCencus,CuratedFusionGene0,[]);

h1=1;
TbyGene_mix_lf = table();
hit_2_include=[];
for c1=1:size(TbyGene_mix,2)
    if (TbyGene_mix(c1).avg_dist == -1 || TbyGene_mix(c1).avg_dist > len_filter)
        hit_2_include(h1)=c1;
        h1=h1+1;
    end
end
TbyGene_mix_lf = TbyGene_mix(hit_2_include);
        

annotated_table = annotate_hits_list( TbyGene_mix_lf,SVTable,bins,hitstable_mix_lookup,pa_mix );
hits_table=table();
hits_table.cluster_num = annotated_table.hit_num;
hits_table.sid = annotated_table.sid;
hits_table.gene_i = annotated_table.gene_i;
hits_table.gene_j = annotated_table.gene_j;
hits_table.subtype = annotated_table.subtype;
hits_table.chr_i = annotated_table.seqnames;
hits_table.pos_i = annotated_table.start;
hits_table.strand_i = annotated_table.strand;
hits_table.chr_j = annotated_table.altchr;
hits_table.pos_j = annotated_table.altpos;
hits_table.strand_j = annotated_table.altstrand;
hits_table.pval = annotated_table.pval;
writetable(hits_table,'sigSV_annot','delimiter','\t')
