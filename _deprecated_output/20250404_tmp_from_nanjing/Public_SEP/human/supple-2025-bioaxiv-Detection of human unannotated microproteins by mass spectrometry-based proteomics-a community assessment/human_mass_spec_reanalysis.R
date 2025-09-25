library("ggplot2")

##### FIGURE SIZE
one.c <- 85 #single column
one.5c <- 144 #1.5 column
two.c <- 174 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9


#proteomapper retrieval functions

get_proteomapper_matches<-function(df_peptides,atlas_data) #identify canonical and noncanonical proteins that a list of peptides maps to according to ProteoMapper
{
	#peptideatlas canononical proteins
	canonical_proteins<-read.csv("input_files/query_guest_20221027-103625.csv",sep=",")

	df_peptides$noncanonical_protein_matches<-""
	df_peptides$canonical_protein_matches<-""

	lpeptides<-gsub("I","L",df_peptides$peptide) #I and L treated the same by proteomapper, so convert all I->L 

	for(i in 1:length(lpeptides))
	{
		all_matches<-atlas_data$protein[which(atlas_data$peptide == lpeptides[i] | atlas_data$peptide == df_peptides$peptide[i])]
		canonical_matches<-intersect(all_matches,canonical_proteins$nextprot_accession)
		if(length(canonical_matches)==0)
		{
			noncanonical_matches<-all_matches
		}
		else
		{
			noncanonical_matches<-all_matches[-which(all_matches%in%canonical_matches)]
		}
		df_peptides[i,]$noncanonical_protein_matches<-paste(noncanonical_matches,sep=",",collapse=",")
		df_peptides[i,]$canonical_protein_matches<-paste(canonical_matches,sep=",",collapse=",")	
	}
	return(df_peptides)
}



############### LISTS of noncanonical PSMs from all studies 

#### Chothani 

#read supplementary data file from Chothani 
mmc9<-read.csv("Chothani/mmc9.csv",header=F)

#get peptide sequences without mos 
unmod_peptides<-gsub('[0-9]+','',gsub('[^[:alnum:]]','',substr(mmc9$V10,3,nchar(mmc9$V10)-2)))

#extract scan numbers from file, which are in different columns for different sets of experiments
scannum<-mmc9$V3
alt_scan<-array()
for(i in 1:length(mmc9[,1]))
{
	if(grepl("scan=",mmc9[i,]$V4))
	{
		alt_scan[i]<-substr(mmc9[i,]$V4,gregexpr('scan=',mmc9[i,]$V4)[[1]][1]+5,10000)
	}
	else
	{
		alt_scan[i]<-(-1)
	}
}

scannum[which(scannum==-1 & alt_scan!=-1)]<-alt_scan[which(scannum==-1 & alt_scan!=-1)]
scannum[1:62]<-mmc9$V4[1:62]

#must remove "TESTPEP" string to get correct filenames 
mgf_files<-gsub('TESTPEP','',mmc9$V1)

df_chothani<-data.frame(
	pmid=35841888,
	curator="Aaron",
	peptides=unmod_peptides,
	orf_name=mmc9$V11,
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="",
	filenames=mgf_files,
	scan_num=scannum,
	charge=0,
	mods="",
	usi=""
)

#set PXDs for the four datasets analyzed in study 
df_chothani$pxd[which(mmc9$V18=="ESC")]<-"PXD006271"
df_chothani$pxd[which(mmc9$V18=="Heart")]<-"PXD006675"
df_chothani$pxd[which(mmc9$V18=="Fib")]<-"PXD003406 to PXD003417"
df_chothani$pxd[which(mmc9$V18=="HUVEC")]<-"PXD003406 to PXD003417"


px_ids<-3406:3417
px_pre<-c("HUVEC_cyt_con","HUVEC_cyt_stim","HUVEC_ne_con","HUVEC_ne_stim","HUVEC_sn_con","HUVEC_sn_stim","NHDF_cyt_con","NHDF_cyt_stim","NHDF_ne_con","NHDF_ne_stim","NHDF_sn_con","NHDF_sn_stim")

for(i in 1:length(px_pre))
{
	df_chothani$pxd[which(substr(df_chothani$filename,1,nchar(px_pre[i]))==px_pre[i])]<-paste("PXD00",px_ids[i],sep="")
}

#3406: HUVEC_cyt_con
#3407: HUVEC_cyt_stim
#3408: HUVEC_ne_con
#3409: HUVEC_ne_stim
#3410: HUVEC_sn_con
#3411: HUVEC_sn_stim
#3412: NHDF_cyt_con
#3413: NHDF_cyt_stim
#3414: NHDF_ne_con
#3415: NHDF_ne_stim
#3416: NHDF_sn_con
#3417: NHDF_sn_stim

#write all peptides to file, to be analyzed with proteomapper
write.table(unique(df_chothani$peptide),"chothani_peptides",row.names=F,quote=F,col.names=F)

#read in proteomapper results file 
chothani_atlas<-read.csv("Chothani/chothani_peptide_atlas.xls",sep="\t")

#for all peptides, determine if they map to canonical or noncanonical proteins in peptideatlas 
df_chothani<-get_proteomapper_matches(df_chothani,chothani_atlas)

#cleaning up data to put in standard format
require("stringr")
df_chothani$scan_num<-str_remove(df_chothani$scan_num,"Locus:1.1.1.")
df_chothani$mods<-substr(mmc9$V10,3,nchar(mmc9$V10)-2)
df_chothani$charge<-mmc9$V9
df_chothani$mods<-str_replace_all(df_chothani$mods, '[0-9+.]+', function(x) paste("[",x,"]",sep=""))
df_chothani$mods<-str_replace_all(df_chothani$mods, coll("[+15.995]"),coll("[OXIDATION]"))
df_chothani$mods<-str_replace_all(df_chothani$mods, coll("C[+57.021]"),coll("C[CARBAMIDOMETHYL]"))
df_chothani$mods<-str_replace_all(df_chothani$mods, coll("[+42.011]"),coll("[ACETYL]"))

#make USIs
df_chothani$usi<-paste("mzspec:",
	df_chothani$pxd,":",
	substr(df_chothani$filenames,1,nchar(df_chothani$filenames)-4),":scan:",df_chothani$scan_num,":",df_chothani$mods,"/",df_chothani$charge,	sep="")


############ Martinez 

#read all PSM files from supplementary data 
X<-vector()
X<-rbind(X,read.csv("Martinez/SupB15-WT_smORF_hits.txt",sep="\t",header=F))
X<-rbind(X,read.csv("Martinez/SupB15-RT_smORF_hits.txt",sep="\t",header=F))
X<-rbind(X,read.csv("Martinez/JY_smORF_hits.txt",sep="\t",header=F))
X<-rbind(X,read.csv("Martinez/HCT116_smORF_hits.txt",sep="\t",header=F))
X<-rbind(X,read.csv("Martinez/HCC1937_smORF_hits.txt",sep="\t",header=F))
X<-rbind(X,read.csv("Martinez/HCC1143_smORF_hits.txt",sep="\t",header=F))
X<-rbind(X,read.csv("Martinez/fibroblast_smORF_hits.txt",sep="\t",header=F))

#process data to get PSM information
all_peptides<-substr(X$V2,3,nchar(X$V2)-2)
all_proteins<-substr(X$V18,2,nchar(X$V18)-2)
protein_info<-data.frame(orf="",peptide=all_peptides,filename=X$V4, scannum=X$V5, charge=X$V3)

reform_coords<-array()
for(i in 1:length(all_proteins))
{
	chrstr<-strsplit(all_proteins[i],":")[[1]][1]
	coords<-strsplit(all_proteins[i],":")[[1]][2]
	coord_list<-strsplit(coords,"-")[[1]]
	reform_coords[i]<-paste("smORF-",chrstr,":",as.numeric(coord_list[1])+2,"-",coord_list[2],sep="")
}


every_protein<-vector()
num_proteins<-array()
protein_info_all<-vector()
for(i in 1:length(all_proteins))
{
	num_proteins[i]<-length(strsplit(all_proteins[i],",")[[1]])
	every_protein<-c(every_protein,strsplit(all_proteins[i],",")[[1]])
	plist<-strsplit(all_proteins[i],",")[[1]]
	for(j in 1:length(plist))
	{
		protein_info_all<-rbind(protein_info_all,protein_info[i,])
		protein_info_all[length(protein_info_all$orf),]$orf<-plist[j]
	}
}



protein_info_all$orf<-trimws(protein_info_all$orf)
protein_info_all<-protein_info_all[which(grepl("smORF",protein_info_all$orf)),]

every_protein<-protein_info_all$orf
all_proteins_mod<-array()
for(i in 1:length(every_protein))
{
	dum<-strsplit(every_protein[i],"-")[[1]]
	if(substr(every_protein[i],6,6)=="+")
	{
		dum1<-as.numeric(strsplit(every_protein[i],"-")[[1]][2])
		all_proteins_mod[i]<-paste(dum[1],"-",dum1+3,sep="")
	}
	else
	{
		dum1<-strsplit(dum[2],":")[[1]]
		all_proteins_mod[i]<-paste("smORF-",dum1[1],":",as.numeric(dum1[2])-2,"-",as.numeric(dum[3]),sep="")
	}
}

protein_info_all$alt_coords<-all_proteins_mod

#read and process second set of PSMs from pFind results 
Y<-vector()
Y<-rbind(Y,read.csv("Martinez/fib_mHLA.pFind.protein_.txt",sep="\t",header=F))
Y<-rbind(Y,read.csv("Martinez/HCC1143.pFind.protein_.txt",sep="\t",header=F))
Y<-rbind(Y,read.csv("Martinez/HCC1937.pFind.protein_.txt",sep="\t",header=F))
Y<-rbind(Y,read.csv("Martinez/HCT116_mHLA.pFind.protein_.txt",sep="\t",header=F))
Y<-rbind(Y,read.csv("Martinez/JY_mHLA.pFind.protein_.txt",sep="\t",header=F))
Y<-rbind(Y,read.csv("Martinez/SupB15_RT_mHLA.pFind.protein_.txt",sep="\t",header=F))
Y<-rbind(Y,read.csv("Martinez/SupB15_WT_mHLA.pFind.protein_.txt",sep="\t",header=F))

sel<-which(grepl("smORF",Y$V11))
pfind_data<-data.frame(proteins=Y$V11[sel], peptide=Y$V4[sel], filename=Y$V25[sel], p1="",p2="",p3="",p4="",p5="",p6="",p7="",p8="",p9="",p10="",charge=Y$V26[sel],mods=Y$V9[sel])

qval0<-as.numeric(Y$V4[sel-1])
miss<-which(is.na(qval0))
i<-1
while(length(miss)>0)
{
	qval0[miss]<-as.numeric(Y$V4[sel[miss]-i-1])
	miss<-which(is.na(qval0))
	i<-i+1
}
pfind_data$qval<-qval0

nummatches<-array()
for(i in 1:length(pfind_data[,1]))
{
	nummatches[i]<-length(strsplit(pfind_data[i,]$proteins,"/")[[1]])
	pfind_data[i,4:(3+nummatches[i])]<-strsplit(pfind_data[i,]$proteins,"/")[[1]]
}

pfind_data$nummatches<-nummatches

#pfind_data_cut<-pfind_data[which(pfind_data$qval<.01),]
pfind_data_cut<-pfind_data[which(pfind_data$qval<1),]

pfind_data_cut_all<-vector()

for(i in 1:length(pfind_data_cut[,1]))
{
	for(j in 1:10)
	{
		if(pfind_data_cut[i,3+j]!="")
		{
			pfind_data_cut_all<-rbind(pfind_data_cut_all,pfind_data_cut[i,])
			pfind_data_cut_all$proteins[length(pfind_data_cut_all$proteins)]<-pfind_data_cut[i,3+j]
		}
	}
}

pfind_data_cut_all<-pfind_data_cut_all[which(grepl("smORF",pfind_data_cut_all$proteins)),]

pfind_data_cut_all$scannum=""
pfind_data_cut_all$file=""

for(i in 1:length(pfind_data_cut_all[,1]))
{
	pfind_data_cut_all[i,]$scannum<-strsplit(pfind_data_cut_all$filename[i],'\\.')[[1]][2]
	pfind_data_cut_all[i,]$file<-strsplit(pfind_data_cut_all$filename[i],'\\.')[[1]][1]
}


martinez_full<-data.frame(orf=c(protein_info_all$orf,pfind_data_cut_all$proteins),peptide=c(protein_info_all$peptide,pfind_data_cut_all$peptide),filename=c(protein_info_all$filename,pfind_data_cut_all$file),scannum=c(protein_info_all$scannum,pfind_data_cut_all$scannum), charge=c(protein_info_all$charge,pfind_data_cut_all$charge),mods=c(rep("",length(protein_info_all$orf)),pfind_data_cut_all$mods))
martinez_full[1021,]$orf<-"smORF+chr1:37940194-37941091"
martinez_full[1068,]$orf<-"smORF-chr4:110223217-110223625"

martinez_full$coord_mod<-array()

detected_list<-read.csv("Martinez/martinez_detected_full.csv",sep=",")
detected_list[157,1]<-"chr1:37940194-37941094"

for(i in 1:length(martinez_full$orf))
{
	dum<-strsplit(martinez_full$orf[i],"-")[[1]]
	if(substr(martinez_full$orf[i],6,6)=="+")
	{
		dum1<-as.numeric(strsplit(martinez_full$orf[i],"-")[[1]][2])
		martinez_full$coord_mod[i]<-paste(dum[1],"-",dum1+3,sep="")
	}
	else
	{
		dum1<-strsplit(dum[2],":")[[1]]
		martinez_full$coord_mod[i]<-paste("smORF-",dum1[1],":",as.numeric(dum1[2])-2,"-",as.numeric(dum[3]),sep="")
	}
}
martinez_full$included<-0
martinez_full$included[which(substr(martinez_full$coord_mod,7,100) %in% detected_list[,1])]<-1


every_protein<-unique(c(pfind_data_cut_all$proteins,protein_info_all$orf))
every_proteins_mod<-array()
for(i in 1:length(every_protein))
{
	dum<-strsplit(every_protein[i],"-")[[1]]
	if(substr(every_protein[i],6,6)=="+")
	{
		dum1<-as.numeric(strsplit(every_protein[i],"-")[[1]][2])
		every_proteins_mod[i]<-paste(dum[1],"-",dum1+3,sep="")
	}
	else
	{
		dum1<-strsplit(dum[2],":")[[1]]
		every_proteins_mod[i]<-paste("smORF-",dum1[1],":",as.numeric(dum1[2])-2,"-",as.numeric(dum[3]),sep="")
	}
}

write.table(unique(martinez_full$peptide),"martinez_peptides",row.names=F,col.names=F,quote=F)

martinez_pmap<-read.csv("Martinez/martinez_proteomapper.txt",sep="\t")

df_martinez<-data.frame(
	pmid=31819274,
	curator="Aaron",
	peptides=martinez_full$peptide,
	orf_name=martinez_full$orf,
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD000394",
	filenames=martinez_full$filename,
	scan_num=martinez_full$scannum,
	charge=martinez_full$charge,
	mods=martinez_full$peptide,
	usi=""
)

df_martinez<-get_proteomapper_matches(df_martinez,martinez_pmap)

for(i in 1:length(martinez_full$mods))
{
	allmods<-rev(strsplit(martinez_full$mods[i],";")[[1]])
	if(length(allmods)>0)
	{
		for(j in 1:length(allmods))
		{
			columns<-strsplit(allmods[j],",")[[1]]
			if(columns[2]=="Gln->pyro-Glu[AnyN-termQ]")
			{
				df_martinez$mods[i]<-paste("[Gln->pyro-Glu]",df_martinez$mods[i],sep="")
			}
			else
			{
				df_martinez$mods[i]<-paste(substr(df_martinez$mods[i],0,as.numeric(columns[1])),"[",substr(columns[2],1,nchar(columns[2])-3),"]",substr(df_martinez$mods[i],as.numeric(columns[1])+1,nchar(df_martinez$mods[i])),sep="")
			}
		}
	}
}

df_martinez<-df_martinez[which(martinez_full$included==1),]

df_martinez$usi<-paste("mzspec:",
	df_martinez$pxd,":",
	df_martinez$filenames,":scan:",
	df_martinez$scan_num,":",
	df_martinez$mods,"/",
	df_martinez$charge,sep="")


########## Ouspenskaia

#read in PSM data
regev_files<-list.files("regev/Tables")
Y<-list()

for(i in 1:length(regev_files))
{
	Y[[i]]<-read.csv(paste("regev/Tables/",regev_files[i],sep=""),sep="\t",header=T)
}

for(i in 1:length(Y))
{
	Y[[i]]<-Y[[i]][which(Y[[i]]$condType=="nuORF"),]
}
all_seq<-vector()
all_orfs<-vector()
all_files<-vector()
num_orfs<-array()
for(i in 1:length(Y))
{
	num_orfs[i]<-length(Y[[i]]$Mapped.protein)
	if(length(Y[[i]]$sequence>0))
	{
		all_seq<-c(all_seq,Y[[i]]$sequence)
		if(i!=28)
		{
			all_orfs<-c(all_orfs,Y[[i]]$Mapped.protein)
		}
		else
		{
			all_orfs<-c(all_orfs,Y[[i]]$accession_numbers_nuORF)
		}
		all_files<-c(all_files,Y[[i]]$filename)
	}
}

all_filenames<-array()
all_scannum<-array()
all_charges<-array()
for(i in 1:length(all_files))
{
	all_filenames[i]<-strsplit(all_files[i],"\\.")[[1]][1]
	all_scannum[i]<-strsplit(all_files[i],"\\.")[[1]][2]
	all_charges[i]<-strsplit(all_files[i],"\\.")[[1]][4]
	
}

df_regev<-data.frame(
	pmid=34663921,
	curator="Aaron",
	peptides=all_seq,
	orf_name=all_orfs,
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="",
	filenames=all_filenames,
	scan_num=all_scannum,
	charge=all_charges,
	mods=all_seq,
	usi=""
)

#sort out which filename is associated with which MSV
df_regev$pxd[which(substr(df_regev$filename,1,5)=="M2015"|substr(df_regev$filename,1,5)=="M2016")]<-"MSV000080527"
df_regev$pxd[which(substr(df_regev$filename,1,5)=="C2015")]<-"MSV000080527"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180507")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,6)=="GN2017")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180516")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180519")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="GG20161008")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180724")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,24)=="YY_20180828_SK_HLA_B5703")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,24)=="YE_20180726_SK_HLA_B1517")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="GG20161104")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180430")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180502")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,6)=="GG2017")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,6)=="AC2017")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="GG20170301")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180511")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180517")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180428")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,3)=="H2A")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,19)=="GG20161007_CRH_HLAC")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YY_20180830")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180504")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180531")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180512")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,6)=="MS2016")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180510")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,2)=="YD")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,5)=="G2014")]<-"MSV000084442"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="YJ20180319")]<-"MSV000084442"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YE_20180622")]<-"MSV000084442"
df_regev$pxd[which(substr(df_regev$filename,1,2)=="FE")]<-"MSV000084442"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="YJ20180308")]<-"MSV000084787"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="GG20160819")]<-"MSV000084787"
df_regev$pxd[which(substr(df_regev$filename,1,22)=="YE_20180726_SK_HLA_Mel")]<-"MSV000084442"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="GG20160930")]<-"MSV000084442"
df_regev$pxd[which(substr(df_regev$filename,1,2)=="YG")]<-"MSV000084787"
df_regev$pxd[which(substr(df_regev$filename,1,6)=="YY2019")]<-"MSV000084787"
df_regev$pxd[which(substr(df_regev$filename,1,17)=="GG20161007_CRH_CP")]<-"Unknown"
df_regev$pxd[which(substr(df_regev$filename,1,10)=="YJ20180316")]<-"MSV000084442"
df_regev$pxd[which(grepl("QE2_B",df_regev$filename))]<-"Unknown"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YY_20180813")]<-"MSV000084172"
df_regev$pxd[which(substr(df_regev$filename,1,11)=="YY_20180817")]<-"MSV000084787"


# regev mods

# MS/MS search parameters included: no-enzyme specificity; fixed modification: cysteinylation of cysteine; variable modifications: carbamidomethylation of cysteine, oxidation of methionine, and pyroglutamic acid at peptide N-terminal glutamine; Variable modification of carbamidomethylation of cysteine was only used for HLA alleles that included an alkylation step

# whole proteome: Carbamidomethylation of cysteine was set as a fixed modification. For the GBM dataset TMT labeling was required at lysine, but peptide N-termini were allowed to be either labeled or unlabeled. Allowed variable modifications were acetylation at the protein N-terminus, oxidized methionine, pyroglutamic acid, deamidated asparagine and pyrocarbamidomethyl cysteine

files0<-list.files("regev/search/PSMresultsExport_B721221_alleles/Alleles_ABC_noIAA",full.names=T)
files1<-list.files("regev/search/PSMresultsExport_B721221_alleles/Alleles_ABCG_IAA",full.names=T)
files2<-list.files("regev/search/PSMresultsExport_patients/PatientsIP_nuORF",full.names=T)
files3<-list.files("regev/search/PSMresultsExport_patients/PatientsIP_snvIndel_nuORF",full.names=T)
files4<-list.files("regev/search/PSMresultsExport_B721221_alleles/Blanks/",full.names=T)

search_files<-c(files0,files1,files2,files3,files4)


regev_search<-read.csv(search_files[1],sep=";")
for(i in 2:length(search_files))
{
	table_read<-read.csv(search_files[i],sep=";")
	if(length(colnames(table_read))==61)
	{
		regev_search<-rbind(regev_search,table_read)
	}
	else
	{
		regev_search<-rbind(regev_search,table_read[,c(1:55,57:62)])
	}
}

#df_regev$mods<-regev_search[file_match,]$sequence
#mods<-regev_search[file_match,]$modification

library("stringr")
df_regev$mods[which(grepl("HLA",df_regev$filenames))]<-str_replace_all(df_regev$mod[which(grepl("HLA",df_regev$filenames))],coll("C"),coll("C[CYSTEINYL]"))
df_regev$mods[which(!grepl("HLA",df_regev$filenames))]<-str_replace_all(df_regev$mod[which(!grepl("HLA",df_regev$filenames))],coll("C"),coll("C[CARBAMIDOMETHYL]"))


df_regev$mods[which(grepl("TMT",df_regev$filenames))]<-str_replace_all(df_regev$mod[which(grepl("TMT",df_regev$filenames))],coll("K"),coll("K[TMT6PLEX]"))
df_regev$mods[which(grepl("TMT",df_regev$filenames))]<-paste("[TMT6PLEX]-",df_regev$mod[which(grepl("TMT",df_regev$filenames))],sep="")

df_regev$mods<-str_replace_all(df_regev$mods,coll("m"),coll("M[OXIDATION]"))
df_regev$mods<-str_replace_all(df_regev$mods,coll("q"),coll("Q[PRO->PYRO-GLU]"))

df_regev<-df_regev[which(all_files  %in% regev_search$filename & df_regev$pxd!="Unknown"),]

df_regev$usi<-paste("mzspec:",
	df_regev$pxd,":",
	df_regev$filename,":scan:",
	df_regev$scan_num,":",
	df_regev$mod,"/",
	df_regev$charge,sep="")

write.table(unique(df_regev$peptide),"regev_peptides",row.names=F,col.names=F,quote=F)



regev_pmap0<-read.csv("regev/regive_atlas0.txt",sep="\t")
regev_pmap1<-read.csv("regev/regive_atlas1.txt",sep="\t")
regev_pmap<-rbind(regev_pmap0,regev_pmap1)

df_regev<-get_proteomapper_matches(df_regev,regev_pmap)

############# Lu

library("stringr")
mq<-read.csv("Lu/Lu_maxquant.csv",skip=2)
xt<-read.csv("Lu/Lu_xtandem.csv",skip=2)
both<-read.csv("Lu/Lu_both.csv",skip=2)

xt_scans<-rep(0,length(xt[,1]))
xt_scans[which(substr(xt$Spectrum.Title,1,12)=="Locus:1.1.1.")]<-substr(xt$Spectrum.Title,13,11+str_locate(substr(xt$Spectrum.Title,13,1000), "\\.")[,1])[which(substr(xt$Spectrum.Title,1,12)=="Locus:1.1.1.")]
xt_scans[which(substr(xt$Spectrum.Title,1,12)!="Locus:1.1.1.")]<-substr(xt$Spectrum,str_locate(xt$Spectrum.Title,"\\.")+1,str_locate(xt$Spectrum.Title,"\\.")[,1]+str_locate(substr(xt$Spectrum,str_locate(xt$Spectrum.Title,"\\.")[,1]+1,1000),"\\.")-1)[which(substr(xt$Spectrum.Title,1,12)!="Locus:1.1.1.")]

both_scans<-rep(0,length(both[,1]))
both_scans[which(substr(both$Spectrum.name,1,12)=="Locus:1.1.1.")]<-substr(both$Spectrum.name,13,11+str_locate(substr(both$Spectrum.name,13,1000), "\\.")[,1])[which(substr(both$Spectrum.name,1,12)=="Locus:1.1.1.")]
both_scans[which(substr(both$Spectrum.name,1,12)!="Locus:1.1.1.")]<-substr(both$Spectrum.name,str_locate(both$Spectrum.name,"\\.")+1,str_locate(both$Spectrum.name,"\\.")[,1]+str_locate(substr(both$Spectrum.name,str_locate(both$Spectrum.name,"\\.")[,1]+1,1000),"\\.")-1)[which(substr(both$Spectrum.name,1,12)!="Locus:1.1.1.")]

xt_filenames<-substr(xt$Spectrum.Title,str_locate(xt$Spectrum.Title,"File")[,1]+6,nchar(xt$Spectrum.Title)-6)
xt_filenames[which(substr(xt$Spectrum.Title,1,12)!="Locus:1.1.1.")]<-substr(xt$Spectrum.Title,1,str_locate(xt$Spectrum.Title,"\\.")-1)[which(substr(xt$Spectrum.Title,1,12)!="Locus:1.1.1.")]

both_filenames<-substr(both$Spectrum.name,str_locate(both$Spectrum.name,"File")[,1]+6,nchar(both$Spectrum.name)-6)
both_filenames[which(substr(both$Spectrum.name,1,12)!="Locus:1.1.1.")]<-substr(both$Spectrum.name,1,str_locate(both$Spectrum.name,"\\.")-1)[which(substr(both$Spectrum.name,1,12)!="Locus:1.1.1.")]


scan_numbers<-c(mq$Scan.number,xt_scans,both_scans)

filenames<-c(mq$Raw.file,xt_filenames,both_filenames)

mq$pxd<-"Unknown"
mq$pxd[which(substr(mq$Raw.file,1,4)=="CHPP")]<-"PXD000529"
mq$pxd[which(substr(mq$Raw.file,1,2)=="3B")]<-"PXD000533"
mq$pxd[which(substr(mq$Raw.file,1,1)=="H" & substr(mq$Raw.file,1,3)!="Hep")]<-"PXD000533"
mq$pxd[which(substr(mq$Raw.file,1,3)=="97H")]<-"PXD000535"
mq$pxd[which(substr(mq$Raw.file,1,3)=="Hep")]<-"PXD000535"
mq$pxd[which(substr(mq$Raw.file,1,3)=="CYZ")]<-"IPX0000762000"
mq$pxd[which(substr(mq$Raw.file,1,4)=="2015")]<-"IPX0000763000"


xt$pxd<-"Unknown"
xt$pxd[which(substr(xt_filenames,1,4)=="CHPP")]<-"PXD000529"
xt$pxd[which(substr(xt_filenames,1,2)=="3B")]<-"PXD000533"
xt$pxd[which(substr(xt_filenames,1,1)=="H" & substr(xt_filenames,1,3)!="Hep")]<-"PXD000533"
xt$pxd[which(substr(xt_filenames,1,3)=="97H")]<-"PXD000535"
xt$pxd[which(substr(xt_filenames,1,3)=="Hep")]<-"PXD000535"
xt$pxd[which(substr(xt_filenames,1,3)=="CYZ")]<-"IPX0000762000"
xt$pxd[which(substr(xt_filenames,1,4)=="2015")]<-"IPX0000763000"

both$pxd<-"Unknown"
both$pxd[which(substr(both_filenames,1,4)=="CHPP")]<-"PXD000529"
both$pxd[which(substr(both_filenames,1,2)=="3B")]<-"PXD000533"
both$pxd[which(substr(both_filenames,1,1)=="H" & substr(both_filenames,1,3)!="Hep")]<-"PXD000533"
both$pxd[which(substr(both_filenames,1,3)=="97H")]<-"PXD000535"
both$pxd[which(substr(both_filenames,1,3)=="Hep")]<-"PXD000535"
both$pxd[which(substr(both_filenames,1,3)=="CYZ")]<-"IPX0000762000"
both$pxd[which(substr(both_filenames,1,4)=="2015")]<-"IPX0000763000"

mq_mods<-substr(mq$Modified.sequence,2,nchar(mq$Modified.sequence)-1)
mq_mods<-str_replace_all(mq_mods,coll("M(ox)"),coll("M[OXIDATION]"))
mq_mods<-str_replace_all(mq_mods,coll("(ac)"),coll("[ACETYL]"))
mq_mods<-str_replace_all(mq_mods,coll("(gl)"),coll("[GLN->PYRO-GLU]"))

mod_reformat<-function(mod_text,mod_symbol,modtable,peptides,shift)
{
	reformatted_peptipdes<-peptides
	oxidation_open<-str_locate(modtable,mod_text)[,1]
	oxidation_close<-str_locate(substr(modtable,oxidation_open,1000),"\\)")[,1]
	#substr(modtable,oxidation_open+nchar(mod_text)+1,oxidation_close-1)
	for(i in 1:length(peptides))
	{
		indices<-as.numeric(rev(strsplit(substr(modtable[i],oxidation_open[i]+nchar(mod_text)+1,oxidation_open[i]+oxidation_close[i]-2),",")[[1]]))
		if(length(indices)>0)
		{
			for(j in 1:length(indices))
			{
				if(!is.na(indices[j]))
				{

					reformatted_peptipdes[i]<-paste(substr(reformatted_peptipdes[i],1,indices[j]+shift),mod_symbol,substr(reformatted_peptipdes[i],indices[j]+1+shift,1000),sep="")
				}
			}
		}
	}
	return(reformatted_peptipdes)
}

xt_mods<-xt$Sequence
xt_mods<-mod_reformat("Oxidation of M","[OXIDATION]",xt$Variable.Modification.s.,xt_mods,0)
xt_mods<-mod_reformat("Acetylation of protein N-term","[ACETYL]",xt$Variable.Modification.s.,xt_mods,-1)
xt_mods<-mod_reformat("Pyrolidone from carbamidomethylated C","[Pyro-carbamidomethyl]",xt$Variable.Modification.s.,xt_mods,-1)
xt_mods<-mod_reformat("Pyrolidone from E","[Glu->pyro-Glu]",xt$Variable.Modification.s.,xt_mods,-1)
xt_mods<-mod_reformat("Pyrolidone from Q","[Gln->pyro-Glu]",xt$Variable.Modification.s.,xt_mods,-1)


both_mods<-both$Peptide.sequence

for(i in 1:length(both_mods))
{
	all_mods<-rev(trimws(strsplit(both$Fixed.modifications.identified.by.spectrum[i],split=",")[[1]]))
	if(length(all_mods)>0)
	{
		for(j in 1:length(all_mods))
		{
			if(grepl("Oxidation",all_mods[j]))
			{
				mod_index<-as.numeric(substr(all_mods[j],2,nchar(all_mods[j])-20))
				both_mods[i]<-paste(substr(both_mods[i],1,mod_index),"[OXIDATION]",substr(both_mods[i],mod_index+1,1000),sep="")
			}
		}
	}
}

for(i in 1:length(both_mods))
{
	all_mods<-rev(trimws(strsplit(both$Fixed.modifications.identified.by.spectrum[i],split=",")[[1]]))
	if(length(all_mods)>0)
	{
		for(j in 1:length(all_mods))
		{
			if(grepl("n-term: Acetyl",all_mods[j]))
			{
				both_mods[i]<-paste("[ACETYL]",both_mods[i],sep="")
			}
			if(grepl("n-term: Gln->pyro-Glu",all_mods[j]))
			{
				both_mods[i]<-paste("[Gln->pyro-Glu]",both_mods[i],sep="")
			}
		}
	}
}

peptides<-toupper(c(mq$Sequence,xt$Sequence,both$Peptide.sequence))
write.table(unique(peptides),"lu_peptides",row.names=F,quote=F,col.names=F)

df_lu<-data.frame(
	pmid=31340039,
	curator="Aaron",
	peptides=peptides,
	orf_name=c(mq$Proteins,xt$Protein.s., both$Protein.name), 
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd=c(mq$pxd,xt$pxd,both$pxd),
	filenames=filenames,
	scan_num=scan_numbers,
	charge=c(mq$Charge,xt$Idenitification.Charge,both$Spectrum.charge),
	mods=c(mq_mods,xt_mods,toupper(both_mods)),
	usi=""
)

df_lu$mods<-str_replace_all(df_lu$mods,coll("[ACETYL]"),coll("[acetyl]"))
df_lu$mods<-str_replace_all(df_lu$mods,coll("C"),coll("C[CARBAMIDOMETHYL]"))
df_lu$mods<-str_replace_all(df_lu$mods,coll("[acetyl]"),coll("[ACETYL]"))

df_lu$usi<-paste("mzspec:",df_lu$pxd,":",df_lu$filenames,":scan:",df_lu$scan_num,":",df_lu$mods,"/",df_lu$charge,sep="")

lu_atlas<-read.csv("Lu/Lu_peptide_atlas.xls.txt",sep="\t")

df_lu<-get_proteomapper_matches(df_lu,lu_atlas)



##### duffy

### from Supplementary Table 2 of Duffy 
prenatal<-read.csv("Duffy/prenatal_brain.csv")
ngn2<-read.csv("Duffy/ngn2_neurons.csv")
adult<-read.csv("Duffy/adult_brain.csv")

noncanonical_types<-c("D-uORF","E-overlap","F-internal","G-external","H-polycistronic","I-readthrough","J-noncoding")

prenatal_peptides<-read.csv("Duffy/peptide_prenatal.txt",sep="\t") # peptides file output by MaxQuant
ngn2_peptides<-read.csv("Duffy/peptides_ngn2_correctversion_try2.txt",sep="\t")
adult_peptides<-read.csv("Duffy/peptides_adult_cortex.txt",sep="\t")

prenatal_evidence<-read.csv("Duffy/evidence_prenatal.txt",sep="\t") 
ngn2_evidence<-read.csv("Duffy/evidence_ngn2.txt",sep="\t") 
adult_evidence<-read.csv("Duffy/evidence_adult_cortex.txt",sep="\t") 

prenatal_filtered<-prenatal_evidence[which(prenatal_evidence$Leading.razor.protein %in% prenatal$orfID[which(prenatal$type %in% noncanonical_types)]),]
ngn2_filtered<-ngn2_evidence[which(ngn2_evidence$Leading.razor.protein %in% ngn2$orfID[which(ngn2$type %in% noncanonical_types)]),]
adult_filtered<-adult_evidence[which(adult_evidence$Leading.razor.protein %in% adult$orfID[which(adult$type %in% noncanonical_types)]),]


peptides<-unique(c(prenatal_filtered$Sequence,ngn2_filtered$Sequence,adult_filtered$Sequence))
write.table(unique(peptides),"duffy_peptides",row.names=F,quote=F,col.names=F)

library("stringr")
mods<-c(prenatal_filtered$Modified.sequence,ngn2_filtered$Modified.sequence,adult_filtered$Modified.sequence)
mods<-substr(mods,2,nchar(mods)-1)
mods<-str_replace_all(mods,coll("C"),coll("C[CARBAMIDOMETHYL]"))
mods<-str_replace_all(mods,coll("M(Oxidation (M))"),coll("M[OXIDATION]"))
mods<-str_replace_all(mods,coll("(Acetyl (Protein N-term))"),coll("[ACETYL]"))



df_duffy<-data.frame(
	pmid=36171426,
	curator="Aaron",
	peptides=c(prenatal_filtered$Sequence,ngn2_filtered$Sequence,adult_filtered$Sequence),
	orf_name=c(prenatal_filtered$Leading.proteins,ngn2_filtered$Leading.proteins, adult_filtered$Leading.proteins), 
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD035950",
	filenames=c(prenatal_filtered$Raw.file, ngn2_filtered$Raw.file, adult_filtered$Raw.file),
	scan_num=c(prenatal_filtered$MS.MS.scan.number,ngn2_filtered$MS.MS.scan.number,adult_filtered$MS.MS.scan.number),
	charge=c(prenatal_filtered$Charge,ngn2_filtered$Charge,adult_filtered$Charge),
	mods=mods,
	usi=""
)


duffy_atlas<-read.csv("Duffy/duffy_atlas.xls.txt",sep="\t")
df_duffy<-get_proteomapper_matches(df_duffy,duffy_atlas)

df_duffy$usi<-paste("mzspec:",df_duffy$pxd,":",df_duffy$filenames,":scan:",df_duffy$scan_num,":",df_duffy$mods,"/",df_duffy$charge,sep="")

##### duffy decoys 

decoys_filtered<-prenatal_evidence[which(grepl("REV",prenatal_evidence$Leading.proteins)),]
write.table(unique(decoys_filtered$Sequence),"duffy_decoys",row.names=F,quote=F,col.names=F)
sel_decoys<-c(1,2,3,4,10,11,12,18,19,22)
decoys_used<-decoys_filtered[sel_decoys,]

library("stringr")
decoy_mods<-c(decoys_used$Modified.sequence)
decoy_mods<-substr(decoy_mods,2,nchar(decoy_mods)-1)
decoy_mods<-str_replace_all(decoy_mods,coll("C"),coll("C[CARBAMIDOMETHYL]"))
decoy_mods<-str_replace_all(decoy_mods,coll("M(Oxidation (M))"),coll("M[OXIDATION]"))

decoy_mods<-str_replace_all(decoy_mods,coll("[CARBAMIDOMETHYL]"),coll("[Carbamidomethyl]"))
decoy_mods<-str_replace_all(decoy_mods,coll("[OXIDATION]"),coll("[Oxidation]"))


df_decoys<-data.frame(
	pmid=36171426,
	curator="Aaron",
	peptides=decoys_used$Sequence,
	orf_name=decoys_used$Leading.proteins, 
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD035950",
	filenames=decoys_used$Raw.file, 
	scan_num=decoys_used$MS.MS.scan.number,
	charge=decoys_used$Charge,
	mods=decoy_mods,
	usi=""
)

df_decoys$usi<-paste("mzspec:",df_decoys$pxd,":",df_decoys$filenames,":scan:",df_decoys$scan_num,":",df_decoys$mods,"/",df_decoys$charge,sep="")
df_decoys$instrument<-"LTQ Orbitrap Velos"
df_decoys$LTQ<-"Yes"

######### Chong 

chong<-read.csv("Chong/chong_hla_psms.csv",skip=1)

write.table(unique(chong$Sequence),"chong_peptides",row.names=F,quote=F,col.names=F)

# No fixed mods. For samples 0D5P, 0NVC and 0MM745, oxidation (M) and phosphorylation (STY) were set as variable modifications; for the remaining samples only oxidation (M) was included as a variable modification.

library("stringr")
mods<-chong$Peptide
mods<-str_replace_all(mods,coll("(Oxidation)"),coll("[OXIDATION]"))
mods<-str_replace_all(mods,coll("(Phospho)"),coll("[PHOSPHO]"))

scan_num<-substr(chong$Spectrum,str_locate(chong$Spectrum,"\\.")+1,str_locate(chong$Spectrum,"\\.")[,1]+str_locate(substr(chong$Spectrum,str_locate(chong$Spectrum,"\\.")[,1]+1,1000),"\\.")-1)
filenames<-substr(chong$Spectrum,1,str_locate(chong$Spectrum,"\\.")-1)

df_chong<-data.frame(
	pmid=32157095,
	curator="Aaron",
	peptides=chong$Sequence,
	orf_name=chong$Gene_name, 
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD013649",
	filenames=filenames,
	scan_num=scan_num,
	charge=chong$Charge,
	mods=mods,
	usi=""
)

df_chong$usi<-paste("mzspec:",df_chong$pxd,":",df_chong$filenames,":scan:",df_chong$scan_num,":",df_chong$mods,"/",df_chong$charge,sep="")

chong_atlas<-read.csv("Chong/chong_peptide_atlas.txt",sep="\t")
df_chong<-get_proteomapper_matches(df_chong,chong_atlas)

#### Chen 

chen_hla<-read.csv("Chen/txt_HLA_sORF/txt_HLA/evidence.txt",sep="\t")
chen_full<-read.csv("Chen/txt_FullProteome_sORF/txt_FullProteome_sORF/evidence.txt",sep="\t")

chen_full_peptides<-read.csv("Chen/txt_FullProteome_sORF/txt_FullProteome_sORF/peptides.txt",sep="\t")


table(chen_full$Leading.razor.protein[which(grepl("new",chen_full$Proteins) & chen_full$Proteins==chen_full$Leading.razor.protein)])
#pg<-read.csv("Chen/txt_FullProteome_sORF/txt_FullProteome_sORF/proteinGroups.txt")

table(chen_full$Proteins[which(grepl("upstream",chen_full$Proteins) & chen_full$Modifications=="Unmodified" )])
table(chen_full$Proteins[which(grepl("new",chen_full$Proteins) & chen_full$Modifications=="Unmodified" & chen_full$Proteins==chen_full$Leading.razor.protein)])
table(chen_full$Proteins[which(grepl("start_overlap",chen_full$Proteins) & chen_full$Modifications=="Unmodified" & chen_full$Proteins==chen_full$Leading.razor.protein)])

which(grepl("RAB34_27044563_143aa",chen_full$Proteins))
which(grepl("MYZAP_57998925_86aa",chen_full$Proteins))

peptide_index<-c(which(grepl("upstream",chen_full$Proteins) & chen_full$Modifications=="Unmodified" ),which(grepl("new",chen_full$Proteins) & chen_full$Modifications=="Unmodified" & chen_full$Proteins==chen_full$Leading.razor.protein),which(grepl("start_overlap",chen_full$Proteins) & chen_full$Modifications=="Unmodified" & chen_full$Proteins==chen_full$Leading.razor.protein),which(grepl("RAB34_27044563_143aa",chen_full$Proteins)),which(grepl("MYZAP_57998925_86aa",chen_full$Proteins)))

library("stringr")

df_chen<-data.frame(
	pmid=32139545,
	curator="Aaron",
	peptides=chen_full$Sequence[peptide_index],
	orf_name=chen_full$Proteins[peptide_index], 
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD014031",
	filenames=chen_full$Raw.file[peptide_index], 
	scan_num=chen_full$MS.MS.scan.number[peptide_index],
	charge=chen_full$Charge[peptide_index],
	mods=chen_full$Sequence[peptide_index],
	usi=""
)

df_chen$usi<-paste("mzspec:",df_chen$pxd,":",df_chen$filenames,":scan:",df_chen$scan_num,":",df_chen$mods,"/",df_chen$charge,sep="")

write.table(unique(df_chen$peptides),"chen_peptides",row.names=F,quote=F,col.names=F)
chen_atlas<-read.csv("Chen/chen_atlas.xls.txt",sep="\t")
df_chen<-get_proteomapper_matches(df_chen,chen_atlas)




### Prensner

library("stringr")
presner<-read.csv("Prensner/presner_peptides.csv",skip=2)
presner<-presner[which(!is.na(presner$scanNumber)),]
presner_filenames<-substr(presner$filename,1,str_locate(presner$filename,"\\.")[,1]-1)

pdc_ids<-read.csv("Prensner/PDC_file_manifest_04212023_180525_.tsv.txt",sep="\t")
pdc_filenames<-substr(pdc_ids$File.Name,1,nchar(pdc_ids$File.Name)-4)

presner_pdc_ids<-pdc_ids$PDC.Study.ID[match(presner_filenames,pdc_filenames)]

mods<-presner$sequence


CPTAC_TMT_SET<-c(33212010,30205044,31031003,31675502,31988290, 32059776,32649874,33212010)


mods[which(presner$Publication..PMID..if.external.dataset.%in% CPTAC_TMT_SET)]<-str_replace_all(mods[which(presner$Publication..PMID..if.external.dataset.%in% CPTAC_TMT_SET)],coll("C"),coll("C[CARBAMIDOMETHYL]"))
mods[which(presner$Publication..PMID..if.external.dataset.%in% c(31844290))]<-str_replace_all(mods[which(presner$Publication..PMID..if.external.dataset.%in% c(31844290))],coll("C"),coll("C[CYSTEINYL]"))

mods[which(grepl("c:carbamidomethylation",presner$varMods))]<-str_replace_all(mods[which(grepl("c:carbamidomethylation",presner$varMods))],coll("c"),coll("C[CARBAMIDOMETHYL]"))
mods[which(grepl("c:pyroCarbamidomethylCys",presner$varMods))]<-str_replace_all(mods[which(grepl("c:pyroCarbamidomethylCys",presner$varMods))],coll("c"),coll("C[PYRO-CARBAMIDOMETHYL]"))

mods<-str_replace_all(mods,coll("m"),coll("M[OXIDATION]"))

mods<-str_replace_all(mods,coll("t"),coll("T[PHOSPHO]"))
mods<-str_replace_all(mods,coll("s"),coll("S[PHOSPHO]"))

mods[which(presner$Publication..PMID..if.external.dataset.%in% CPTAC_TMT_SET)]<-str_replace_all(mods[which(presner$Publication..PMID..if.external.dataset.%in% CPTAC_TMT_SET)],coll("K"),coll("K[TMT6PLEX]"))

mods<-str_replace_all(mods,coll("k"),coll("K[ACETYL]"))
mods<-str_replace_all(mods,coll("p"),coll("P[OXIDATION]"))
mods<-str_replace_all(mods,coll("n"),coll("N[DEAMIDATED]"))
mods<-str_replace_all(mods,coll("q"),coll("Q[PRO->PYRO-GLU]"))

mods[which(substr(mods,1,1)!="[" & presner$Publication..PMID..if.external.dataset. %in% CPTAC_TMT_SET)]<-paste("[TMT6PLEX]-",mods[which(substr(mods,1,1)!="[" & presner$Publication..PMID..if.external.dataset%in% CPTAC_TMT_SET)],sep="")

#79.98

df_presner<-data.frame(
	pmid=33510483,
	curator="Aaron",
	peptides=toupper(presner$sequence),
	orf_name=presner$ORF.Name,
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="NA",
	filenames=presner_filenames,
	scan_num=presner$scanNumber,
	charge=presner$parent_charge,
	mods=mods,
	usi=""
	)
	
df_presner$pxd[which(presner$Publication==30205044)]<-"MSV000082644"
df_presner$pxd[which(presner$Publication==31844290)]<-"MSV000084172"
df_presner$pxd[which(presner$dataset=="HLA_Ouspenskaia")]<-"MSV000084787"

df_presner$pxd[which(presner$Data.deposition=="https://cptac-data-portal.georgetown.edu/study-summary/S045")]<-"PXD999945" # don't have on peptideatlas
df_presner$pxd[which(presner$Data.deposition=="https://cptac-data-portal.georgetown.edu/study-summary/S050")]<-"PXD999950" # don't have on peptideatlas
df_presner$pxd[which(presner$Data.deposition=="https://cptac-data-portal.georgetown.edu/study-summary/S051")]<-"PXD999951"
df_presner$pxd[which(presner$Data.deposition=="https://cptac-data-portal.georgetown.edu/study-summary/S053")]<-"PXD999953"
df_presner$pxd[which(presner$Data.deposition=="https://cptac-data-portal.georgetown.edu/study-summary/S056")]<-"PXD999956" #don't have on peptideatlas
df_presner$pxd[which(presner$Data.deposition=="https://cptac-data-portal.georgetown.edu/study-summary/S060")]<-"PXD999960"

df_presner$usi<-paste("mzspec:",df_presner$pxd,":",df_presner$filenames,":scan:",df_presner$scan_num,":",df_presner$mods,"/",df_presner$charge,sep="")

write.table(unique(df_presner$peptides),"presner_peptides",row.names=F,quote=F,col.names=F)

presner_atlas<-read.csv("Prensner/presner_atlas.txt",sep="\t")
df_presner<-get_proteomapper_matches(df_presner,presner_atlas)

### Bogaert results 

bog<-read.csv("Bogaert/Bogaert_results.csv",skip=1)
bog<-bog[which(!is.na(bog$Start)),]


df_bogr<-data.frame(
	pmid=35788065,
	curator="Aaron",
	peptides=bog$Sequence,
	orf_name=bog$Accession,
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD030601",
	filenames=substr(bog$Filename,1,str_locate(bog$Filename,"_b")[,1]-1),
	scan_num="",
	charge=bog$Charge,
	mods=bog$Modified.Sequence,
	usi=""
	)


for(i in 1:length(df_bogr$peptides))
{
	df_bogr$scan_num[i]<-strsplit(bog$Filename[i],"_")[[1]][11]
}

#Mods 
#heavy acetylation of lysine side chains (with 13C1D3-acetate), carbamidomethylation of cysteine and methionine oxidation to methionine sulfoxide were set as fixed modifications. Variable modifications were acetylation of N termini (both light and heavy because of the 13C1D3 label) and pyroglutamate formation of Nt-glutamine (both at the peptide level). 

library("stringr")
df_bogr$mods<-str_replace_all(df_bogr$mods, coll("<Mox*>"),coll("[OXIDATION]"))
df_bogr$mods<-str_replace_all(df_bogr$mods, coll("<AcD4*>"),coll("[+46.0298]")) #should be heavy
df_bogr$mods<-str_replace_all(df_bogr$mods, coll("AcD4-"),coll("[+46.0298]-")) #should be heavy
df_bogr$mods<-str_replace_all(df_bogr$mods, coll("<Cmm*>"),coll("[CARBAMIDOMETHYL]")) 
df_bogr$mods<-str_replace_all(df_bogr$mods, coll("Ace-"),coll("[ACETYL]-")) #should be light
df_bogr$mods<-str_replace_all(df_bogr$mods, coll("-COOH"),coll("")) #should be light


df_bogr$usi<-paste("mzspec:",df_bogr$pxd,":",df_bogr$filenames,":scan:",df_bogr$scan_num,":",df_bogr$mods,"/",df_bogr$charge,sep="")

write.table(unique(df_bogr$peptides),"bogr_peptides",row.names=F,quote=F,col.names=F)

bog_atlas<-read.csv("Bogaert/bog_atlas.xls.txt",sep="\t")

df_bogr<-get_proteomapper_matches(df_bogr,bog_atlas)

bogaert_passed_peptides<-c("MDGEEKTCGGCEGPDAMYVKLISSDGHEFIVKR","MKEETKEDAEEKQ","ADDAGAAGGPGGPGGPEMGNRGGFRGGF","FPSNWNEIVDSFDDMNLSESLLR","MASAASSSSLE","GDVVPKDANAAIATIKTKR","PKDANAAIATIKTKR")
df_bogr<-df_bogr[which(df_bogr$peptides %in% bogaert_passed_peptides),]

#### Cao results

cao<-read.csv("Cao/cao_peptides_fixed.csv") #table given to me by paper authors. there are 28 peptides listed here, but according to table S1 in the paper, some of these are n-terminal extensions rather than independent noncanonical ORFs and so we need to exlude those 
peptide_s1_table<-read.csv("Cao/cao_s1_ms_peptides.csv",sep=",") #Supplementary proteomics data table S1 from paper 

n_terminal_extension_peptides<-c("PVSLLAPLIPPR","AAPELGPGATIEAGAAR","AAPELGPGATIEAGAAR","AAPELGPGATIEAGAAR","ADAGAMAATDIAR","AAAAEEAAAAGPR","AAAAAAAAAAAAAAGTR")

df_cao<-data.frame(
	pmid=35393574,
	curator="Aaron",
	peptides=cao$peptide,
	orf_name=peptide_s1_table$Peptide.ID[match(cao$peptide,peptide_s1_table$Peptide.detected.by.shotgun.proteomics)],
	noncanonical_protein_matches="",
	canonical_protein_matches="",
	pxd="PXD026880",
	filenames=cao$file,
	scan_num=cao$scan,
	charge=2,
	mods=cao$peptide,
	usi=""
	)




df_cao$mods[which(grepl("Oxidation \\(M\\)",cao$modification))]<-str_replace_all(df_cao$mods[which(grepl("Oxidation \\(M\\)",cao$modification))],coll("M"),coll("M[OXIDATION]"))

df_cao$mods[which(grepl("Acetyl \\(N-term\\)",cao$modification))]<-paste("[ACETYL]-",df_cao$mods[which(grepl("Acetyl \\(N-term\\)",cao$modification))],sep="")
df_cao$mods[2]<-"[ACETYL]-AAPELGPGATIEAGAAR"

df_cao$charge<-round(cao$Mr.expt./cao$Observed)

df_cao$usi<-paste("mzspec:",df_cao$pxd,":",df_cao$filenames,":scan:",df_cao$scan_num,":",df_cao$mods,"/",df_cao$charge,sep="")
#df_cao<-df_cao[-which(df_cao$peptides %in% n_terminal_extension_peptides),]

write.table(unique(df_cao$peptides),"cao_peptides",row.names=F,quote=F,col.names=F)

cao_atlas<-read.csv("Cao/cao_atlas.xls.txt",sep="\t")
df_cao<-get_proteomapper_matches(df_cao,cao_atlas)
df_cao$orf_name<-paste("cao_orf_",1:length(df_cao$orf_name),sep="")


###Douka
valid_smorfs<-read.csv("input_files/20221005_validation_smORF.csv")###data gathered from supplementary tables 

df_douka<-data.frame(
	pmid=34193551,
	curator="Jana",
	peptides=valid_smorfs[which(valid_smorfs$PMID==34193551),]$Peptide.Sequence,
	orf_name=valid_smorfs[which(valid_smorfs$PMID==34193551),]$ORF.name,
	noncanonical_protein_matches=valid_smorfs[which(valid_smorfs$PMID==34193551),]$reported.as.non.canonical..novel.isoform..pseudo.gene..microprotein.,
	canonical_protein_matches=valid_smorfs[which(valid_smorfs$PMID==34193551),]$Maps.to.Ordinary.protein..enter.ordinary.protein.ID..comma.delimited.,
	pxd=valid_smorfs[which(valid_smorfs$PMID==34193551),]$ProteomeXchange.ID,
	filenames=valid_smorfs[which(valid_smorfs$PMID==34193551),]$ProteomeXchange.ID,
	scan_num=valid_smorfs[which(valid_smorfs$PMID==34193551),]$Filename,
	charge=valid_smorfs[which(valid_smorfs$PMID==34193551),]$Parent.charge,
	mods=valid_smorfs[which(valid_smorfs$PMID==34193551),]$Peptide.with.modification,
	usi=valid_smorfs[which(valid_smorfs$PMID==34193551),]$USI..
	)

df_douka$charge<-c(3,3,3,2,2,2,3,2,2,3,2,2,3,2,2,2,2,2) #17th and 18th not present
df_douka$usi<-paste(df_douka$usi,df_douka$charge,sep="")
df_douka$usi[12]<-"mzspec:PXD014381:OE05115:scan:9669:TC[CARBAMIDOMETHYL]GNKASAGHEATEK/2"

##van Heesch
df_heesch<-data.frame(
	pmid=31155234,
	curator="Jana",
	peptides=valid_smorfs[which(valid_smorfs$PMID==31155234),]$Peptide.Sequence,
	orf_name=valid_smorfs[which(valid_smorfs$PMID==31155234),]$ORF.name,
	noncanonical_protein_matches=valid_smorfs[which(valid_smorfs$PMID==31155234),]$reported.as.non.canonical..novel.isoform..pseudo.gene..microprotein.,
	canonical_protein_matches=valid_smorfs[which(valid_smorfs$PMID==31155234),]$Maps.to.Ordinary.protein..enter.ordinary.protein.ID..comma.delimited.,
	pxd=valid_smorfs[which(valid_smorfs$PMID==31155234),]$ProteomeXchange.ID,
	filenames=valid_smorfs[which(valid_smorfs$PMID==31155234),]$Filename,
	scan_num=valid_smorfs[which(valid_smorfs$PMID==31155234),]$Scan.num,
	charge=valid_smorfs[which(valid_smorfs$PMID==31155234),]$Parent.charge,
	mods=valid_smorfs[which(valid_smorfs$PMID==31155234),]$Peptide.with.modification,
	usi=valid_smorfs[which(valid_smorfs$PMID==31155234),]$USI..
	)
df_heesch<-df_heesch[which(df_heesch$scan_num!=""),]
df_heesch$mods<-str_replace_all(df_heesch$mods,coll("C"),coll("C[CARBAMIDOMETHYL]"))

df_heesch$usi<-paste("mzspec:",df_heesch$pxd,":",df_heesch$filenames,":scan:",df_heesch$scan_num,":",df_heesch$mods,"/",df_heesch$charge,sep="")



#### combining PSM data from all studies 
df_all<-rbind(df_chothani,df_martinez, df_regev, df_lu, df_duffy, df_chong, df_chen, df_presner, df_bogr, df_cao, df_douka, df_heesch)

##we can only consider PSMs from repositories with PXD or MSV ids  
df_all<-df_all[which(substr(df_all$pxd,1,3) %in% c("PXD","MSV")),]
require("stringr")

#standardize capitalization of modifications
df_all$usi<-str_replace_all(df_all$usi, coll("[CARBAMIDOMETHYL]"),coll("[Carbamidomethyl]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[OXIDATION]"),coll("[Oxidation]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[PHOSPHO]"),coll("[Phospho]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[ACETYL]"),coll("[Acetyl]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[TMT6PLEX]"),coll("[TMT6plex]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[CYSTEINYL]"),coll("[Cysteinyl]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[PRO->PYRO-GLU]"),coll("[Pro->pyro-Glu]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[DEAMIDATED]"),coll("[Deamidated]"))
df_all$usi<-str_replace_all(df_all$usi, coll("[Carbamidomethyl]-"),coll("[Carbamidomethyl]"))

## set HLAs, mostly by PMIDs since all studies but Prensner are only all HLA or no HLA; otherwise by filename 
df_all$HLA<-FALSE
df_all$HLA[which(grepl("HLA",df_all$filenames) & df_all$pmid==34663921)]<-TRUE
df_all$HLA[which(grepl("HLA",df_all$filenames) & df_all$pmid==33510483)]<-TRUE

df_all$HLA[which(df_all$pmid==32157095)]<-TRUE
df_all$HLA[which(df_all$pmid==31819274)]<-TRUE

#### PXD TO INSTRUMENT # read in instrument information for each study and at it to PSM dataframe 

instrument_table<-read.table("input_files/instrument_table.txt",sep=";")
df_all$instrument<-""
df_all$LTQ<-"No"

for(i in 1:length(instrument_table[,1]))
{
	df_all$instrument[which(df_all$pxd==instrument_table[i,1])]<-instrument_table[i,2]
	df_all$LTQ[which(df_all$pxd==instrument_table[i,1])]<-instrument_table[i,3]
}

write.csv(df_all,"all_peptide_usis_")


### ploting Figure 1A : orf counts in study vs. reported detections 

studies<-read.csv("input_files/studies_table_large.csv")

df_orfs_vs_detects<-data.frame(num_orfs=as.numeric(studies$Total.number.ORFs),num_detects=studies$Reported.noncanonical.ORFs.with.MS.support,database_source=studies$ORF.database.source,cite=paste(studies$Citation,studies$Year),pmid=studies$PMID)
df_orfs_vs_detects<-df_orfs_vs_detects[order(df_orfs_vs_detects$pmid),]
df_orfs_vs_detects$cite<-factor(df_orfs_vs_detects$cite,levels=df_orfs_vs_detects$cite)

df_orfs_vs_detects$num_orfs[which(df_orfs_vs_detects$database_source=="three frame translation")]<-1000000
df_orfs_vs_detects$author<-substr(df_orfs_vs_detects$cite,0,nchar(as.character(df_orfs_vs_detects$cite))-5)
df_orfs_vs_detects$y_offset<-c(160,100,100,5,140,50,4,-1000,7,3,200,130)

require("ggplot2")
fig_orfs_vs_findings_outlegend<-ggplot(df_orfs_vs_detects,aes(x=num_orfs,y=num_detects,col=database_source))+
	geom_point()+
	theme_bw()+
	scale_x_continuous(trans='log10',breaks=c(100,1000,10000,100000),limits=c(10,460000), labels=c("100\n","1000\n","10,000\n","100,000\n"))+
	scale_y_continuous(trans='log10',limits=c(1,5000))+
	scale_color_manual(values=c("black","#1A85FF","#D41159","orange"))+
	geom_text(aes(label=author,y=num_detects+y_offset),size=2,col="black")+
	labs(x="sORFs in database",y="sORF-encoded protein detections",col="sORF database source")+
	theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size = txt),
	legend.text = element_text(size = txt*.75),
	legend.title = element_text(size=titles*.75),
	plot.margin = unit(c(1,1,1,1), "pt"),
	legend.key.size = unit(.2, "cm"),
	legend.box.background=element_rect(),
	legend.position=c(.2,.85))
	
fig_orfs_vs_findings_outlegend2<-ggplot(df_orfs_vs_detects[df_orfs_vs_detects$num_orfs==1000000,],aes(x=num_orfs,y=num_detects,col=database_source))+
	geom_point()+
	theme_bw()+
	scale_x_continuous(trans='log10',breaks=c(1000000),labels=c("Full\ntranscriptome"))+
	scale_y_continuous(trans='log10',limits=c(1,5000))+
	scale_color_manual(values=c("orange","#1A85FF","#D41159","black"))+
	geom_text(aes(label=author,y=num_detects+y_offset),size=2,col="black")+
	labs(x="",y="",col="Database source")+
	theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size = txt),
	legend.text = element_text(size = txt*.75),
	legend.title = element_text(size=titles*.75),
	plot.margin = unit(c(1,1,1,-1), "pt"),
	legend.key.size = unit(.2, "cm"),
	legend.position="none",
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())
	
require("cowplot")
fig_orfs_vs_findings_combo<-plot_grid(fig_orfs_vs_findings_outlegend,fig_orfs_vs_findings_outlegend2,rel_widths=c(.8,.2),ncol=2,nrow=1)

png("fig_orfs_vs_findings_combo.png",height=75,width=180,res=300,unit="mm")
fig_orfs_vs_findings_combo
dev.off()


### plotting Figure 1C: canonical mapping via proteomapper 

#################################################
## Getting canonincal and noncanonical matches 
################################################

get_proteo_matches<-function(df_peptides,atlas_data,max_diffs)
{
	canonical_proteins<-read.csv("/home/acwach/HumanMS/query_guest_20221027-103625.csv",sep=",")
	peptide_diffs<-array()
	for(i in 1:length(atlas_data[,1]))
	{
		peptide_diffs[i]<-length(which(strsplit(atlas_data[i,]$peptide,"")[[1]]!=strsplit(atlas_data[i,]$in_fasta,"")[[1]]))
	}



	atlas_data_filtered<-atlas_data[which(peptide_diffs<=max_diffs),]
	noncanonical_protein_matches<-array()
	canonical_protein_matches<-array()
	tryptic_matches<-array()

	lpeptides<-gsub("I","L",df_peptides$peptide)

	final_nuc<-substr(atlas_data_filtered$in_fasta,nchar(atlas_data_filtered$in_fasta),nchar(atlas_data_filtered$in_fasta))
	tryptic_peptides<-which(final_nuc %in% c("R","K") & atlas_data_filtered$prevAA %in% c("R","K"))

	for(i in 1:length(lpeptides))
	{
		all_matches<-atlas_data_filtered$protein[which(atlas_data_filtered$peptide == lpeptides[i] | atlas_data_filtered$peptide == df_peptides$peptide[i])]
		canonical_matches<-intersect(all_matches,canonical_proteins$nextprot_accession)
		canonical_match_ids<-which((atlas_data_filtered$peptide == lpeptides[i] | atlas_data_filtered$peptide == df_peptides$peptide[i]) & atlas_data_filtered$protein %in% canonical_proteins$nextprot_accession)
		if(length(canonical_matches)==0)
		{
			noncanonical_matches<-all_matches
		}
		else
		{
			noncanonical_matches<-all_matches[-which(all_matches%in%canonical_matches)]
		}
		noncanonical_protein_matches[i]<-paste(noncanonical_matches,sep=",",collapse=",")
		canonical_protein_matches[i]<-paste(canonical_matches,sep=",",collapse=",")	
		tryptic_match_ids<-intersect(canonical_match_ids,tryptic_peptides)
		if(length(tryptic_match_ids)>0) #(i %in% tryptic_peptides)
		{
			tryptic_matches[i]<-paste(canonical_matches,sep=",",collapse=",")	
		}
		else
		{
			tryptic_matches[i]<-""
		}
	}
	df_results<-data.frame(noncanonical_matches=noncanonical_protein_matches,canonical_matches=canonical_protein_matches,tryptic_matches=tryptic_matches)
	return(df_results)
}

get_canmatch_counts<-function(df_study,study_atlas)
{


	matches0<-get_proteo_matches(df_study,study_atlas,0)
	matches1<-get_proteo_matches(df_study,study_atlas,1)
	matches100<-get_proteo_matches(df_study,study_atlas,100)

	canonical_matches2016<-array()
	canonical_matches2023<-array()
	tryptic_matches<-array()

	canonical_matches2016_1<-array()
	canonical_matches2023_1<-array()
	tryptic_matches_1<-array()

	canonical_matches2016_100<-array()
	canonical_matches2023_100<-array()
	tryptic_matches_100<-array()

	for(i in 1:length(matches0$canonical_matches))
	{
		canonical_matches2016[i]<-length(which(strsplit(matches0$canonical_matches[i],",")[[1]] %in% uniprot2016[,1]))
		canonical_matches2023[i]<-length(strsplit(matches0$canonical_matches[i],",")[[1]])
		tryptic_matches[i]<-length(which(strsplit(matches0$tryptic_matches[i],",")[[1]] %in% uniprot2016[,1]))

		canonical_matches2016_1[i]<-length(which(strsplit(matches1$canonical_matches[i],",")[[1]] %in% uniprot2016[,1]))
		canonical_matches2023_1[i]<-length(strsplit(matches1$canonical_matches[i],",")[[1]])
		tryptic_matches_1[i]<-length(which(strsplit(matches1$tryptic_matches[i],",")[[1]] %in% uniprot2016[,1]))

		canonical_matches2016_100[i]<-length(which(strsplit(matches100$canonical_matches[i],",")[[1]] %in% uniprot2016[,1]))
		canonical_matches2023_100[i]<-length(strsplit(matches100$canonical_matches[i],",")[[1]])
		tryptic_matches_100[i]<-length(which(strsplit(matches100$tryptic_matches[i],",")[[1]] %in% uniprot2016[,1]))
	}

	can_proteins<-unique(strsplit(paste(matches1$canonical_matches[which(canonical_matches2016_1==0)],sep=",",collapse=","),",")[[1]])
	#can_proteins<-unique(strsplit(paste(df_study$canonical_protein_matches[which(canonical_matches2016_100==0 & df_study$canonical_protein_matches!="")],sep=",",collapse=","),",")[[1]])
	can_proteins<-can_proteins[can_proteins!=""]
	
	recently_canon<-can_proteins[which(!(can_proteins %in% uniprot2016[,1]))]
	
	# total ORFs 
	# ORFs supported only by peptides that map to canonical ORFs as of 2023
	# ORFs supported only by peptides that map to canonical ORFs as of 2016
	# ORFs supported only by peptides with proteomapper matches to 2016 canonical ORFs with 1 or fewer SNPs
	# ORFs supported only by peptides with proteomapper matches to 2016 canonical ORFs with 0 SNPs
	# ORFs supported only by tryptic peptides with proteomapper matches to 2016 canonical ORFs with 0 SNPs

	# total PSMs 
	# PSMs that map to canonical ORFs as of 2023
	# PSMs that map to canonical ORFs as of 2016
	# PSMs with proteomapper matches to 2016 canonical ORFs with 1 or fewer SNPs
	# PSMs with proteomapper matches to 2016 canonical ORFs with 0 SNPs
	# PSMs involving tryptic peptides with proteomapper matches to 2016 canonical ORFs with 0 SNPs

	total_orfs<-length(unique(df_study$orf_name))
	can_orfs_2023<-total_orfs-length(unique(df_study$orf_name[which(df_study$canonical_protein_matches=="")]))
	can_orfs_2023_1<-total_orfs-length(unique(df_study$orf_name[which(canonical_matches2023_1==0)]))

	can_orfs_2016<-total_orfs-length(unique(df_study$orf_name[which(canonical_matches2016_100==0)]))
	can_orfs_2016_1<-total_orfs-length(unique(df_study$orf_name[which(canonical_matches2016_1==0)]))
	can_orfs_2016_0<-total_orfs-length(unique(df_study$orf_name[which(canonical_matches2016==0)]))
	can_orfs_tryp<-total_orfs-length(unique(df_study$orf_name[which(tryptic_matches==0)]))

	total_psms<-length(df_study$orf_name)
	can_psms_2023_1<-length(which(canonical_matches2023_1>0))
	can_psms_2023<-length(which(canonical_matches2023_100>0))
	can_psms_2016<-length(which(canonical_matches2016_100>0))
	can_psms_2016_1<-length(which(canonical_matches2016_1>0))
	can_psms_2016_0<-length(which(canonical_matches2016>0))
	can_psms_tryp<-length(which(tryptic_matches>0))

	total_peptides<-length(unique(df_study$peptides))
	can_peptides_2023<-length(unique(df_study$peptides[which(df_study$canonical_protein_matches!="")]))
	can_peptides_2023_1<-length(unique(df_study$peptides[which(canonical_matches2023_1>0)]))
	
	can_peptides_2016<-length(unique(df_study$peptides[which(canonical_matches2016_100>0)]))
	can_peptides_2016_1<-length(unique(df_study$peptides[which(canonical_matches2016_1>0)]))
	can_peptides_2016_0<-length(unique(df_study$peptides[which(canonical_matches2016>0)]))
	can_peptides_tryp<-length(unique(df_study$peptides[which(tryptic_matches>0)]))
	can_peptides_tryp1<-length(unique(df_study$peptides[which(tryptic_matches_1>0)]))

	df_can<-data.frame(total_orfs=total_orfs,can_orfs_2023=can_orfs_2023,can_orfs_2023_1=can_orfs_2023_1,can_orfs_2016=can_orfs_2016,can_orfs_2016_1=can_orfs_2016_1,can_orfs_2016_0=can_orfs_2016_0,can_orfs_tryp=can_orfs_tryp,total_psms=total_psms,can_psms_2023=can_psms_2023,can_psms_2023_1=can_psms_2023_1,can_psms_2016=can_psms_2016,can_psms_2016_1=can_psms_2016_1,can_psms_2016_0=can_psms_2016_0,can_psms_tryp=can_psms_tryp,
	total_peptides=total_peptides,can_peptides_2023=can_peptides_2023,can_peptides_2023_1=can_peptides_2023_1 ,can_peptides_2016=can_peptides_2016,can_peptides_2016_1=can_peptides_2016_1,can_peptides_2016_0=can_peptides_2016_0,can_peptides_tryp=can_peptides_tryp,can_peptides_tryp1=can_peptides_tryp1)

	df_can$can_peptides_nontryp_0<-df_can$can_peptides_2016_0-df_can$can_peptides_tryp
	df_can$can_peptides_nontryp_1<-df_can$can_peptides_2016_1-df_can$can_peptides_tryp1
	df_can$can_peptides_tryp_exact1<-df_can$can_peptides_tryp1-df_can$can_peptides_tryp
	df_can$can_peptides_nontryp_exact1<-df_can$can_peptides_nontryp_1-df_can$can_peptides_nontryp_0
	df_can$can_peptides100<-df_can$can_peptides_2016-df_can$can_peptides_2016_1
	df_can$no_match <- df_can$total_peptides-df_can$can_peptides_tryp-df_can$can_peptides_tryp_exact1-df_can$can_peptides_nontryp_0-df_can$can_peptides_nontryp_exact1


	return(list(df_can,recently_canon))
}

uniprot2016<-read.csv("input_files/human_prots_2016")

douka_atlas<-read.csv("Douka/douka_atlas.xls.txt",sep="\t")
heesch_atlas<-read.csv("/home/acwach/HumanMS/heesch_atlas.xls.txt",sep="\t")

cao_canmatch<-get_canmatch_counts(df_cao,cao_atlas)  #note: don't know orf identify of these peptides
bog_canmatch<-get_canmatch_counts(df_bogr,bog_atlas)
chothani_canmatch<-get_canmatch_counts(df_chothani,chothani_atlas)
duffy_canmatch<-get_canmatch_counts(df_duffy,duffy_atlas)
#cai_canmatch<-get_canmatch_counts(df_cai,cai_atlas)
douka_canmatch<-get_canmatch_counts(df_douka,douka_atlas)

presner_canmatch_nonhla<-get_canmatch_counts(df_all[which(df_all$HLA==FALSE & df_all$pmid==33510483),],presner_atlas)
presner_canmatch_hla<-get_canmatch_counts(df_presner[which(df_all$HLA==TRUE & df_all$pmid==33510483),],presner_atlas)
presner_canmatch<-get_canmatch_counts(df_presner,presner_atlas)

regev_canmatch<-get_canmatch_counts(df_regev,regev_pmap)
chen_canmatch<-get_canmatch_counts(df_chen,chen_atlas)
chong_canmatch<-get_canmatch_counts(df_chong,chong_atlas)
martinez_canmatch<-get_canmatch_counts(df_martinez,martinez_pmap)
heesch_canmatch<-get_canmatch_counts(df_heesch,heesch_atlas)
lu_canmatch<-get_canmatch_counts(df_lu,lu_atlas)

all_canmatch<-rbind(cao_canmatch[[1]],bog_canmatch[[1]],chothani_canmatch[[1]],duffy_canmatch[[1]],douka_canmatch[[1]],presner_canmatch[[1]],regev_canmatch[[1]],chen_canmatch[[1]],chong_canmatch[[1]],martinez_canmatch[[1]],heesch_canmatch[[1]],lu_canmatch[[1]])

all_canmatch$hla_peptides_0<-0
all_canmatch$hla_peptides_0[c(7,9,10)]<-(all_canmatch$can_peptides_nontryp_0+all_canmatch$can_peptides_tryp)[c(7,9,10)]
all_canmatch$hla_peptides_1<-0
all_canmatch$hla_peptides_1[c(7,9,10)]<-(all_canmatch$can_peptides_tryp_exact1+all_canmatch$can_peptides_nontryp_exact1)[c(7,9,10)]
all_canmatch$can_peptides_nontryp_0[c(7,9,10)]<-0
all_canmatch$can_peptides_tryp[c(7,9,10)]<-0
all_canmatch$can_peptides_tryp_exact1[c(7,9,10)]<-0
all_canmatch$can_peptides_nontryp_exact1[c(7,9,10)]<-0

all_canmatch$hla_peptides_0[6]<-(presner_canmatch_hla[[1]]$can_peptides_nontryp_0+presner_canmatch_hla[[1]]$can_peptides_tryp)
all_canmatch$hla_peptides_1[6]<-(presner_canmatch_hla[[1]]$can_peptides_tryp_exact1+presner_canmatch_hla[[1]]$can_peptides_nontryp_exact1)
all_canmatch$can_peptides_nontryp_0[6]<-presner_canmatch_nonhla[[1]]$can_peptides_nontryp_0
all_canmatch$can_peptides_tryp[6]<-presner_canmatch_nonhla[[1]]$can_peptides_tryp
all_canmatch$can_peptides_tryp_exact1[6]<-presner_canmatch_nonhla[[1]]$can_peptides_tryp_exact1
all_canmatch$can_peptides_nontryp_exact1[6]<-presner_canmatch_nonhla[[1]]$can_peptides_nontryp_exact1
all_canmatch$total_peptides[6]<-presner_canmatch_hla[[1]]$total_peptides+presner_canmatch_nonhla[[1]]$total_peptides
all_canmatch$no_match[6]<-presner_canmatch_hla[[1]]$no_match+presner_canmatch_nonhla[[1]]$no_match

rownames(all_canmatch)<-c("Cao 2022","Bogaert 2022","Chothani 2022","Duffy 2022","Douka 2021","Prensner 2021","Ouspenskaia 2021","Chen 2020","Chong 2020","Martinez 2020","van Heesch 2019","Lu 2019") 

df_canmatch<-data.frame(publication=rownames(all_canmatch),total_peptides=all_canmatch$total_peptides,matched_peptides=c(all_canmatch$can_peptides_tryp,all_canmatch$can_peptides_tryp_exact1,all_canmatch$can_peptides_nontryp_0,all_canmatch$can_peptides_nontryp_exact1,all_canmatch$no_match,all_canmatch$hla_peptides_0,all_canmatch$hla_peptides_1),match_type=c(rep("Tryptic perfect match",12),rep("Tryptic 1 SNP",12),rep("Nontryptic perfect match",12),rep("Nontryptic 1 SNP",12),rep("No match",12),rep("HLA perfect match",12),rep("HLA 1 SNP",12) ))
df_canmatch$match_type<-factor(df_canmatch$match_type,levels=c("No match","Tryptic perfect match","Tryptic 1 SNP","Nontryptic perfect match","Nontryptic 1 SNP","HLA perfect match","HLA 1 SNP"))
df_canmatch$match_type<-factor(df_canmatch$match_type,levels=c("Tryptic perfect match","Tryptic 1 SNP","Nontryptic perfect match","Nontryptic 1 SNP","HLA perfect match","HLA 1 SNP","No match"))

df_canmatch$prop<-df_canmatch$matched_peptides/df_canmatch$total_peptides

publist_labels<-df_orfs_vs_detects$cite
df_canmatch$publication<-factor(df_canmatch$publication,levels=rownames(all_canmatch))

df_canmatch$publication<-factor(df_canmatch$publication,levels=publist_labels)

fig_canmatch<-ggplot(df_canmatch,aes(x=publication,y=100*matched_peptides/total_peptides, fill=match_type))+
	geom_bar(stat="identity", position="stack")+#position_dodge2())+
	theme_bw()+
	#geom_rect(ymin=1.01,ymax=1.06,xmin=.55,xmax=4.5,col="black",fill="white")+
	#scale_y_continuous(limits=c(0,1.03))+#,breaks=1:5)+
	#scale_fill_manual(values=c("azure4","blue","cyan","darkgreen","chartreuse","brown","red"))+
	scale_fill_manual(values=c("azure4","black","darkgreen","chartreuse","brown","red","royalblue2"))+
	#geom_text(x = 2.5, y = 1.035, label = "HLA",size=3.5)+
	#geom_rect(ymin=1.01,ymax=1.06,xmin=4.5,xmax=13.5,col="black",fill="white")+
	#geom_text(x = 9, y = 1.035, label = "Non-HLA",size=3.5)+
	
	#annotate("text", x = 10, y = 5.36, label = "non-HLA",size=3.5)+
	scale_x_discrete(guide = guide_axis(angle=45),labels=publist_labels)+
	labs(x="",y="Percent peptides mapping\nto canonical proteins")+
	theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size = txt),
	plot.margin = unit(c(10,1,1,10), "pt"),
	legend.title=element_blank(),
	legend.text = element_text(size = txt*.75),
	legend.key.size = unit(.2, "cm"),
	legend.box.margin=margin(-10,0,-10,-10))

	#theme(legend.title=element_blank(),plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))
png("fig_canmatch_.png",height=150,width=150,res=300,unit="mm")
plot(fig_canmatch)
dev.off()


#### Adding columns to supplementary table giving canonical matches 

df_all$has_canonical_match_peptideatlas<-FALSE
df_all$has_canonical_match_uniprotkb2016<-FALSE
df_all$has_canonical_match_peptideatlas[which(df_all$canonical_protein_matches!="")]<-TRUE

uniproktb2016_canonical_matches<-array()
for(i in 1:length(df_all[,1]))
{
	uniproktb2016_canonical_matches[i]<-length(which(strsplit(df_all$canonical_protein_matches[i],",")[[1]] %in% uniprot2016[,1]))
}
df_all$has_canonical_match_uniprotkb2016[which(uniproktb2016_canonical_matches>0)]<-TRUE

### plotting Figure 1B/D: overlap between studies in reported detections 

unique_peptides<-unique(df_all$peptides)
studies_per_peptide<-array()
for(i in 1:length(unique_peptides))
{
	studies_per_peptide[i]<-length(unique(df_all$pmid[which(df_all$peptides==unique_peptides[i])]))
}

# of 9414 peptides, only 326 found in more than one study
table(studies_per_peptide)
length(studies_per_peptide)
length(which(studies_per_peptide>1))/length(studies_per_peptide) #3.5% found more than once 


pmids<-unique(df_all$pmid)
pmids<-sort(pmids)
singleton_fraction<-array()
for(i in 1:length(pmids))
{
	pmid_peptides<-unique(df_all$peptides[which(df_all$pmid==pmids[i])])
	singleton_fraction[i]<-length(which(studies_per_peptide[which(unique_peptides %in% pmid_peptides)]==1))/length(which(studies_per_peptide[which(unique_peptides %in% pmid_peptides)]>=1))
}

df_singleton<-data.frame(cite=df_orfs_vs_detects$cite, multiple_fraction=1-singleton_fraction)
library("ggplot2")
f1d<-ggplot(df_singleton,aes(x=cite,y=multiple_fraction*100))+
	geom_bar(stat="identity",fill="black")+
	theme_bw()+
	scale_x_discrete(guide = guide_axis(angle=45))+
	#labs(x="",y="Percent reported noncanonical\npeptides found in another study")+
	labs(x="",y="Percent reported sORF-encoded\npeptides found in another study")+
	theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size = txt),
	plot.margin = unit(c(10,1,1,10), "pt"))

	#geom_text(aes(label=study_count,x=study_pmid,y=study_count+30))
	
png("f1d.png",height=125,width=90,res=300,unit="mm")
plot(f1d)
dev.off()

unique_peptides_unmapped<-unique(df_all$peptides[which(df_all$canonical_protein_matches=="")])
studies_per_peptide_unmapped<-array()
for(i in 1:length(unique_peptides_unmapped))
{
	studies_per_peptide_unmapped[i]<-length(unique(df_all$pmid[which(df_all$peptides==unique_peptides_unmapped[i])]))
}

pmids<-unique(df_all$pmid)
pmids<-sort(pmids)
singleton_fraction_unmapped<-array()
for(i in 1:length(pmids))
{
	pmid_peptides<-unique(df_all$peptides[which(df_all$pmid==pmids[i] & df_all$canonical_protein_matches=="")])
	singleton_fraction_unmapped[i]<-length(which(studies_per_peptide_unmapped[which(unique_peptides_unmapped %in% pmid_peptides)]==1))/length(which(studies_per_peptide_unmapped[which(unique_peptides_unmapped %in% pmid_peptides)]>=1))
}

df_singleton_unmapped<-data.frame(cite=df_orfs_vs_detects$cite, multiple_fraction=1-singleton_fraction_unmapped)
library("ggplot2")
f1e<-ggplot(df_singleton_unmapped,aes(x=cite,y=multiple_fraction*100))+
	geom_bar(stat="identity",fill="black")+
	theme_bw()+
	scale_x_discrete(guide = guide_axis(angle=45))+
	labs(x="",y="Percent reported sORF-encoded\npeptides found in another study\nexcluding ProteoMapper hits")+
	theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size = txt),
	plot.margin = unit(c(10,1,1,10), "pt"))

	#geom_text(aes(label=study_count,x=study_pmid,y=study_count+30))
	
png("f1e.png",height=125,width=90,res=300,unit="mm")
plot(f1e)
dev.off()




library("cowplot")
Fig1ac<-plot_grid(fig_orfs_vs_findings_combo,fig_canmatch,ncol=1,nrow=2,labels=c('A','C'),label_size=10)#,align="h")
Fig1bd<-plot_grid(f1d,f1e,ncol=1,nrow=2,labels=c('B','D'),label_size=10,align="h")

Fig1<-plot_grid(Fig1ac,Fig1bd,ncol=2,nrow=1,rel_widths=c(1.8,1))

png("Figure1_noa_.png",width=174,height=144,res=300,units = "mm")
Fig1
dev.off()



### Figure 2a: judge agreement 

############## analysis of expert opinion

expert1<-read.csv("filled_lists/PSM__list1.csv")
expert1<-expert1[which(expert1$Instrument!=""),]
expert1$curator<-"A"

expert2<-read.csv("filled_lists/PSM__list2.csv")
expert2<-expert2[which(expert2$Instrument!=""),]
expert2$curator<-"B"

expert3<-read.csv("filled_lists/PSM__list3.csv")
expert3<-expert3[which(expert3$Instrument!=""),]
expert3$curator<-"C"

expert4<-read.csv("filled_lists/PSM__list4.csv")
expert4<-expert4[which(expert4$Instrument!=""),]
expert4$curator<-"D"

expert5<-read.csv("filled_lists/PSM__list5.csv")
expert5<-expert5[which(expert5$Instrument!=""),]
expert5$curator<-"E"

expert6<-read.csv("filled_lists/PSM__list6.csv")
expert6<-expert6[which(expert6$Instrument!=""),]
expert6$curator<-"F"

experts<-rbind(expert5[,-6],expert2[,-6],expert3,expert6[,-c(6:11)],expert4[,-6],expert1[,-c(6,7)])

experts_match<-match(experts$USI,df_all$usi)
experts$pmid<-df_all$pmid[experts_match]
experts$HLA<-df_all$HLA[experts_match]
experts$pmid[which(experts$USI %in% df_decoys$usi)]<-"Control"
experts$Rating..1.5.<-as.numeric(experts$Rating..1.5.)
experts$peptide<-df_all$peptides[experts_match]
experts$orf<-df_all$orf_name[experts_match]

experts$pmid_hla<-experts$pmid
experts$pmid_hla[which(experts$pmid==33510483 & experts$HLA==TRUE)]<-"33510483_HLA"
experts<-experts[-which(is.na(experts$pmid)),]

df_all$assigned<-FALSE
df_all$assigned[which(df_all$usi %in% experts$USI)]<-TRUE

#### importantt numbers 

length(unique(experts$USI[experts$pmid!="Control" & !is.na(experts$Rating)])) # total evaluated USIs (406)
length(unique(experts$peptide[experts$pmid!="Control" & !is.na(experts$Rating)])) # total evaluated peptides (307)
length(unique(experts$orf[experts$pmid!="Control" & !is.na(experts$Rating)])) # total evaluated ORFs/proteins (204)




doubles<-names(which(table(experts$USI)==2))

r1<-array()
r2<-array()
for(i in 1:length(doubles))
{
	sel<-which(experts$USI == doubles[i])

	r1[i]<-as.numeric(experts[sel,]$Rating[1])
	r2[i]<-as.numeric(experts[sel,]$Rating[2])
}
isna<-which(is.na(r1) | is.na(r2))
r1<-r1[-isna]
r2<-r2[-isna]



df_doubles<-data.frame(judge1=integer(),judge2=integer(),hitcount=integer())
for(i in 1:5)
{
	for(j in i:5)
	{
		df_doubles<-rbind(df_doubles,c(judge1=i,judge2=j,hitcount=length(which( (r1==i & r2==j)|(r1==j & r2==i) ))))
	}
}
colnames(df_doubles)<-c("judge1","judge2","hitcount")



p2<-ggplot(df_doubles,aes(x=judge1,y=judge2,size=hitcount))+
	geom_point(aes(size=hitcount))+
	scale_size(range=c(2,9)/2)+
	annotate("text", x = 3.5, y = 2, label = "r=0.82",size=2.5)+
	#scale_shape_manual(values=c(15, 20))+
	#scale_color_manual(values=c("blue","orange"))+
	theme_bw()+
	labs(x="Lower rating",y="Higher rating",size="Count")+
	theme(plot.margin = unit(c(5.5,-5.5,5.5,5.5), "pt"),axis.title.x=element_text(vjust=0,size=titles),axis.title.y=element_text(size=titles),
	axis.text = element_text(size = txt),
	legend.text = element_text(size = txt*.75),
	legend.title = element_text(size=titles*.75),
	legend.box.margin=margin(0,0,0,-10))

### Figure 2B: study ratings 


curators<-unique(experts$curator)
pmids_hla<-unique(experts$pmid_hla)
rating_means<-array(dim=c(length(pmids_hla),length(curators)))
ratings_se<-array(dim=c(length(pmids_hla),length(curators)))

df_all$pmid_hla<-df_all$pmid 
df_all$pmid_hla[which(df_all$pmid==33510483 & df_all$HLA==TRUE)]<-"33510483_HLA"
for(i in 1:length(pmids_hla))
{
	for(j in 1:length(curators))
	{
		rating_means[i,j]<-mean(experts$Rating..1.5.[which(experts$pmid_hla==pmids_hla[i] & experts$curator==curators[j])],na.rm=T)
		n<-length(which(experts$pmid_hla==pmids_hla[i] & experts$curator==curators[j] & !is.na(experts$Rating..1.5.)))
		N<-length(which(df_all$pmid_hla==pmids_hla[i] & df_all$canonical_protein_matches=="")) 
		finite_pop_modification<-sqrt((N-n)/(N-1))
		ratings_se[i,j]<-finite_pop_modification*sd(experts$Rating..1.5.[which(experts$pmid_hla==pmids_hla[i] & experts$curator==curators[j])],na.rm=T)/sqrt(n)
	}
}

publist<-c("van Heesch 2019","Douka 2021","Duffy 2022","Negative control","Bogaert 2022","Ouspenskaia 2021","Cao 2022","Martinez 2020","Chong 2020","Chothani 2022","Prensner 2021 non-HLA","Lu 2019","Prensner 2021 HLA","Chen 2020")
df_study_ratings<-data.frame(ratings=as.vector(rating_means),ratings_se=as.vector(ratings_se), publication=publist,judge=c(rep("A",14),rep("B",14),rep("C",14),rep("D",14),rep("E",14),rep("F",14) ))

publist_labels<-c(as.character(df_orfs_vs_detects$cite),"Negative control")
publist_labels<-c(publist_labels[1:5],"Prensner 2021 non-HLA","Prensner 2021 HLA",publist_labels[7:13])
df_study_ratings$publication<-factor(df_study_ratings$publication,levels=publist_labels)

library(viridis)

fig_study_ratings<-ggplot(df_study_ratings[-which(is.na(df_study_ratings$ratings_se)),],aes(x=publication,y=ratings))+
	geom_errorbar(aes(x=publication,ymin=ratings-ratings_se,ymax=ratings+ratings_se,col=judge),position=position_dodge(.6),width=.3,linewidth = .5) +
	geom_point(aes(x=publication,y=ratings,col=judge),position=position_dodge(.6),size=1.4)+ # ,position=position_dodge(.6),width=.3,linewidth = 1.3) +
	
	geom_bar(stat="identity", position="dodge")+#position_dodge2())+
	theme_bw()+
	scale_color_viridis(discrete = TRUE)+ #scale_color_manual(values = color_values)+
	#geom_rect(ymin=5.2,ymax=5.5,xmin=.45,xmax=4.5,col="black",fill="white")+
	#annotate("text", x = 3, y = 5.36, label = "HLA",size=3.5)+ 
	#geom_rect(ymin=5.2,ymax=5.5,xmin=4.5,xmax=14.5,col="black",fill="white")+
	#annotate("text", x = 10, y = 5.36, label = "non-HLA",size=3.5)+ 
	geom_hline(yintercept=4, linetype='dotted')+
	scale_x_discrete(guide = guide_axis(angle=45),labels=publist_labels )+
	scale_y_continuous(limits=c(.75,5.4),breaks=1:5)+
	labs(x="",y="Mean PSM rating",color="Evaluator",fill="")+
	#theme(plot.margin = unit(c(1,1,1,10), "pt"),
	theme(plot.margin = unit(c(1,1,-10,10), "pt"),
	axis.title.x=element_text(size=titles),
	axis.title.y=element_text(size=titles),
	axis.text = element_text(size = txt),
	legend.text = element_text(size = txt*.75),
	legend.title = element_text(size=titles),
	legend.box.margin=margin(0,0,0,-10)
	)

### Figure 2C: rating distribution 
se_pro<-function(p,n){return(sqrt(p*(1-p)/n))}

HLA_ratings_freqs<-table(experts$Rating..1.5.[experts$pmid!="Control" & experts$HLA==TRUE])
HLA_ratings_freqs<-HLA_ratings_freqs/sum(HLA_ratings_freqs)
non_HLA_ratings_freqs<-table(experts$Rating..1.5.[experts$pmid!="Control" & experts$HLA==FALSE])
non_HLA_ratings_freqs<-non_HLA_ratings_freqs/sum(non_HLA_ratings_freqs)

N_hla<-length(which(experts$HLA==TRUE & !is.na(experts$Rating..1.5.)))
N_nonhla<-length(which(experts$HLA==FALSE & !is.na(experts$Rating..1.5.)))

stderr_hla<-array()
stderr_nonhla<-array()
for(i in 1:5)
{
	stderr_hla[i]<-se_pro(HLA_ratings_freqs[i],N_hla)
	stderr_nonhla[i]<-se_pro(HLA_ratings_freqs[i],N_nonhla)

}
df_rating_hit<-data.frame(scores=c(1:5,1:5),rates=c(HLA_ratings_freqs,non_HLA_ratings_freqs),rates_se=c(stderr_hla,stderr_nonhla), is_hla<-c(rep("HLA",5),rep("non-HLA",5)))


p3<-ggplot(df_rating_hit,aes(x=scores,y=rates,fill=is_hla))+
	geom_bar(stat="identity",position="dodge")+
	geom_errorbar(position="dodge",aes(x=scores,ymin=rates-rates_se,ymax=rates+rates_se)) +
	scale_fill_manual(values=c("#7E909A","#AC3E31"))+
	theme_bw()+
	theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt"),
	axis.title.x=element_text(vjust=0,size=titles),
	axis.title.y=element_text(size=titles),
	axis.text=element_text(size=txt),legend.key.size = unit(.2, "cm"),
	legend.position=c(.45,.8),
	legend.text = element_text(size = txt*.85),
	legend.title = element_text(size=titles*.75),
	)+
	labs(x="PSM score",y="Frequency across all\nstudies/evaluators",fill="")


### Figure 4d/e: ORF analysis: ribo-seq reads and length  

cand<-read.csv("input_files/candidate_orfs",sep=" ")
tcalls<-read.csv("input_files/translation_calls",sep=" ")

special_ids<-unique(cand$special_id)[2:length(unique(cand$special_id))]
read_sum<-array()
special_rating<-array()
special_orflength<-array()
special_pvals<-array()
for(i in 1:length(special_ids))
{
	sel<-which(cand$special_id==special_ids[i])[1]
	read_sum[i]<-tcalls$reads0[sel]
	special_rating[i]<-mean(experts[which(experts$USI==special_ids[i]),]$Rating)
	special_orflength[i]<-cand$ORF_length[sel]
	special_pvals[i]<-1-pbinom(tcalls[sel,]$frame0,tcalls[sel,]$frame_sum,1/3)
}

read_sum<-read_sum[which(!is.na(special_rating))]
special_orflength<-special_orflength[which(!is.na(special_rating))]
special_pvals<-special_pvals[which(!is.na(special_rating))]
special_rating<-special_rating[which(!is.na(special_rating))]


df_special<-data.frame(spec_rating=special_rating,spec_orflength=special_orflength,spec_reads=log(read_sum+.0001),spec_readrate=log((read_sum+.0001)/(special_orflength/3) ))
df_special<-df_special[which(df_special$spec_rating!=3),]
df_special$rated<-"High\nrated"
df_special$rated[which(df_special$spec_rating<3)]<-"Low\nrated"


png("special_orflength.png",width=144,height=144,res=300,units = "mm")
fig_spec_orflength<-ggplot(df_special,aes(x=rated,y=spec_orflength/3-1))+
geom_boxplot(outlier.shape = NA)+
geom_jitter(color="black", size=0.4, alpha=0.9)+
geom_signif(comparisons = list(c("High\nrated", "Low\nrated")), 
              annotations="p=.01",textsize=3, vjust=-1, y_position=170,tip_length=.01)+
coord_cartesian(ylim=c(0,240))+
labs(x="",y="Protein length (aa)")+
theme_bw()+
theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size = txt),
	plot.margin = unit(c(0,5.5,-10,5.5), "pt"))
fig_spec_orflength
dev.off()

library("ggsignif")
png("special_reads.png",width=144,height=144,res=300,units = "mm")
fig_spec_reads<-ggplot(df_special,aes(x=rated,y=spec_readrate))+
geom_boxplot(outlier.shape = NA)+
geom_jitter(color="black", size=0.4, alpha=0.9)+
geom_signif(comparisons = list(c("High\nrated", "Low\nrated")), 
              annotations="p=.005",textsize=3,vjust=-1)+
coord_cartesian(ylim=c(0, 11))+
labs(x="",y="Expression level\n(log in-frame read count per codon)")+
theme_bw()+
theme(plot.title = element_blank(),
	axis.title = element_text(size = titles),
	axis.text = element_text(size =txt))
fig_spec_reads
dev.off()


### figure 4f: extrapolation 

pmids_hla<-unique(experts$pmid_hla)
ratings_hi_count<-array()
for(i in 1:length(pmids_hla))
{
	ratings_hi_count[i]<-length(unique(experts$orf[which(experts$pmid_hla==pmids_hla[i] & experts$Rating>3)]))
}

publist_hla<-c("van Heesch 2019","Douka 2021","Duffy 2022","Negative control","Bogaert 2022","Ouspenskaia 2021 (HLA)","Cao 2022","Martinez 2020 (HLA)","Chong 2020 (HLA)","Chothani 2022","Prensner 2021","Lu 2019","Presner 2021 (HLA)","Chen 2020")

total_orfs_pmid<-array()
tested_orfs_pmid<-array()
for(i in 1:length(pmids_hla))
{
	total_orfs_pmid[i]<-length(unique(df_all$orf_name[(which(df_all$pmid==pmids_hla[i]))]))
	tested_orfs_pmid[i]<-length(unique(df_all$orf_name[(which(df_all$pmid==pmids_hla[i] & df_all$assigned==TRUE))]))
}

total_orfs_pmid[11]<-length(unique(df_all$orf_name[(which(df_all$pmid==33510483 & df_all$HLA==FALSE))]))
tested_orfs_pmid[11]<-length(unique(df_all$orf_name[(which(df_all$pmid==33510483 & df_all$HLA==FALSE & df_all$assigned==TRUE))]))

total_orfs_pmid[13]<-length(unique(df_all$orf_name[(which(df_all$pmid==33510483 & df_all$HLA==TRUE))]))
tested_orfs_pmid[13]<-length(unique(df_all$orf_name[(which(df_all$pmid==33510483 & df_all$HLA==TRUE & df_all$assigned==TRUE))]))

total_orfs_pmid[7]<-17
tested_orfs_pmid[7]<-17

df_hicount<-data.frame(hi_count=ratings_hi_count, publication=publist_hla,HLA=grepl("HLA",publist_hla),total_orfs=total_orfs_pmid,tested_orfs=tested_orfs_pmid)
df_hicount$HLA[which(df_hicount$HLA==TRUE)]<-"HLA"
df_hicount$HLA[which(df_hicount$HLA==FALSE)]<-"non-HLA"
df_hicount$HLA[which(df_hicount$publication=="Negative control")]<-"control"

df_hicount$predicted_orfs<-(df_hicount$hi_count/df_hicount$tested_orfs)*df_hicount$total_orfs

df_hicount$publication<-factor(df_hicount$publication,levels=c("Ouspenskaia 2021 (HLA)","Martinez 2020 (HLA)","Chong 2020 (HLA)","Prensner 2021 (HLA)","van Heesch 2019","Douka 2021","Duffy 2022","Bogaert 2022","Cao 2022","Chothani 2022","Prensner 2021","Lu 2019","Chen 2020","Negative control"))

df_hicount$publication<-factor(df_hicount$publication,levels=c("Chong 2020 (HLA)","Martinez 2020 (HLA)","Ouspenskaia 2021 (HLA)","Prensner 2021 (HLA)","Bogaert 2022","Cao 2022","Chen 2020","Chothani 2022","Douka 2021","Duffy 2022","Lu 2019","Prensner 2021","van Heesch 2019","Negative control"))


#high rated orfs in nonhla studies
hi_rated_nonhla<-sum(df_hicount$hi_count[which(df_hicount$HLA=="non-HLA")])

#tested orfs in nonhla studies 
tested_nonhla<-sum(df_hicount$tested_orfs[which(df_hicount$HLA=="non-HLA")])

#extrapolated high rated orfs in nonhla studies
extrap_nonhla<-sum(df_hicount$predicted_orfs[which(df_hicount$HLA=="non-HLA")])

#total orfs in nonhla studies
total_orfs_nonhla<-sum(df_hicount$total_orfs[which(df_hicount$HLA=="non-HLA")])

#high rated orfs in hla studies
hi_rated_hla<-sum(df_hicount$hi_count[which(df_hicount$HLA=="HLA")])

#tested orfs in hla studies 
tested_hla<-sum(df_hicount$tested_orfs[which(df_hicount$HLA=="HLA")])

#extrapolated high rated orfs in nonhla studies
extrap_hla<-round(sum(df_hicount$predicted_orfs[which(df_hicount$HLA=="HLA")]))

#total orfs in nonhla studies
total_orfs_hla<-sum(df_hicount$total_orfs[which(df_hicount$HLA=="HLA")])

df_nonhla_hicount<-data.frame(orf_counts=c(hi_rated_nonhla,tested_nonhla-hi_rated_nonhla),orf_class=c("High rated","Low rated"))

df_nonhla_extrap<-data.frame(orf_counts=c(round(extrap_nonhla),round(total_orfs_nonhla-extrap_nonhla)),orf_class=c("High rated","Low rated"))

df_hla_hicount<-data.frame(orf_counts=c(hi_rated_hla,tested_hla-hi_rated_hla),orf_class=c("High rated","Low rated"))

df_hla_extrap<-data.frame(orf_counts=c(extrap_hla,total_orfs_hla-extrap_hla),orf_class=c("High rated","Low rated"))

df_nonhla_tested<-rbind(df_nonhla_hicount,df_nonhla_extrap)
df_nonhla_tested$status<-c("Evaluated","Evaluated","Extrapolated","Extrapolated")
df_nonhla_tested[2,1]<-sum(df_nonhla_hicount$orf_counts)
df_nonhla_tested[4,1]<-sum(df_nonhla_extrap$orf_counts)
df_nonhla_tested$orf_class[which(df_nonhla_tested$orf_class=="Low rated")]<-"All"
df_nonhla_tested$rating<-df_nonhla_tested$orf_class
df_nonhla_tested$hla<-"Non-HLA"
df_nonhla_tested$orf_sd<-0
df_nonhla_tested$orf_sd[3]<-27.5

df_hla_tested<-rbind(df_hla_hicount,df_hla_extrap)
df_hla_tested$status<-c("Evaluated","Evaluated","Extrapolated","Extrapolated")
df_hla_tested[2,1]<-sum(df_hla_hicount$orf_counts)
df_hla_tested[4,1]<-sum(df_hla_extrap$orf_counts)
df_hla_tested$orf_class[which(df_hla_tested$orf_class=="Low rated")]<-"All (HLA)"
df_hla_tested$orf_class[which(df_hla_tested$orf_class=="High rated")]<-"High rated (HLA)"

df_hla_tested$rating<-c("High rated","All","High rated","All")
df_hla_tested$hla<-"HLA"


df_hla_tested$orf_sd<-0
df_hla_tested$orf_sd[3]<-29

df_tested<-rbind(df_nonhla_tested,df_hla_tested)

fig_line_extrap<-ggplot(df_tested[which(df_tested$rating!="All"),],aes(x=status,y=orf_counts,group=orf_class,col=hla))+
	geom_line()+#aes(linetype=rating))+
	geom_point()+
	scale_y_continuous(trans='log10')+
	scale_alpha_continuous(range = c(0, 1),guide = 'none')+
	theme_bw()+
	scale_fill_manual(values=c("#7E909A","#AC3E31"))+
	labs(y="High-rated protein detections",x="")+
	geom_errorbar(aes(ymin=orf_counts-orf_sd,ymax=orf_counts+orf_sd,alpha=orf_sd/20),width=.4) +
	theme(legend.title = element_blank(),
	axis.title = element_text(size=titles),
	axis.text = element_text(size = txt),
	legend.text = element_text(size = txt*.75),
	legend.position=c(.25,.8),
	plot.margin = unit(c(0,5.5,-10,5.5),"pt")
		)
	
png("line_extrap.png",width=120,height=80,res=300,unit="mm")
plot(fig_line_extrap)
dev.off()
	


top_fig2<-plot_grid(p2,fig_study_ratings,ncol=2,nrow=1,label_size=10,labels=c('A','B'),rel_widths=c(2,3))
bot_fig2<-plot_grid(p3,fig_spec_reads,fig_spec_orflength,fig_line_extrap,ncol=4,nrow=1,label_size=10,labels=c('C','D','E','F'),rel_widths=c(.75,.5,.5,1))

png("Figure2_updated2.png",width=174,height=144,res=300,units = "mm")
plot_grid(top_fig2,bot_fig2,ncol=1,nrow=2)
dev.off()


############## plots for machine learning spectra comparison 


library(ggplot2)
library(tidyr)

# 1. Create Stacked Bar Plot of Rank vs SA

data <- read.table("Rank_SA_Counts", header = TRUE, sep = "\t")

long_data <- pivot_longer(
  data,
  cols = -Rank,
  names_to = "Value",
  values_to = "Count"
)

long_data$Value <- factor(
  long_data$Value,
  levels = c("High", "Moderate", "Poor", "Terrible", "No.Match")
)

custom_colors <- c(
  "R1" = "red",
  "R2" = "orange",
  "R3" = "yellow",
  "R4" = "lightgreen",
  "R5" = "darkgreen"
)

# Create the bar plot
fig_ml_plot1<-ggplot(long_data, aes(x = Value, y = Count, fill = Rank)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Curated Rating vs Similarity to Predicted Spectra",
    x = "Spectral Similarity to Predicted Spectra",
    y = "Number of Spectra",
    fill = "Curated Rating"
  )

png("m1_plot1.png")
fig_ml_plot1
dev.off()

# 2. Create Box Plot of Rank vs SA

data <- read.table("Rank_SA_Values", header = TRUE, sep = "\t")

fig_ml_plot2<-ggplot(data, aes(x = Curated.Rating, y = Best.Spectral.Angle, fill = Curated.Rating)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2) +
  theme_bw() +
  scale_fill_manual(values = c("R1" = "red", "R2" = "orange", "R3" = "yellow", "R4" = "lightgreen", "R5" = "darkgreen")) +
  scale_x_discrete(labels=1:5)+
  labs(
    y = "Best Spectral Angle",
    x = "Evaluator rating"
  ) + 
  theme(
	legend.position="none",
	axis.text = element_text(size = txt),
	axis.title=element_text(size=titles),
	plot.margin = unit(c(5.5,5.5,5.5,15.5),"pt")
	)

png("m1_plot2.png")
fig_ml_plot2
dev.off()


library("cowplot")
new_top_fig2<-plot_grid(p2,fig_ml_plot2,ncol=2,nrow=1,label_size=10,labels=c('A','B'),rel_widths=c(1,1))
row3<-plot_grid(p3,fig_spec_reads,ncol=2,nrow=1,label_size=10,labels=c('D','E'),rel_widths=c(1.5,1))
row4<-plot_grid(fig_spec_orflength,fig_line_extrap,ncol=2,nrow=1,label_size=10,labels=c('F','G'),rel_widths=c(1,1.5))

png("Figure2_new.png",width=174,height=230,res=300,units = "mm")
plot_grid(new_top_fig2,fig_study_ratings,row3,row4,ncol=1,nrow=4,label_size=10,labels=c("","C",""),rel_heights=c(1,1.5,1,1))
dev.off()



# top_fig2<-plot_grid(p2,fig_study_ratings,ncol=2,nrow=1,label_size=10,labels=c('A','B'),rel_widths=c(2,3))
# bot_fig2<-plot_grid(p3,fig_spec_reads,fig_spec_orflength,fig_line_extrap,ncol=4,nrow=1,label_size=10,labels=c('C','D','E','F'),rel_widths=c(.75,.5,.5,1))

# png("Figure2_updated2.png",width=174,height=144,res=300,units = "mm")
# plot_grid(top_fig2,bot_fig2,ncol=1,nrow=2)
# dev.off()



# p2<-ggplot(df_doubles,aes(x=judge1,y=judge2,size=hitcount))+
	# geom_point(aes(size=hitcount))+
	# scale_size(range=c(2,9)/2)+
	# annotate("text", x = 3.5, y = 2, label = "r=0.82",size=2.5)+
	# #scale_shape_manual(values=c(15, 20))+
	# #scale_color_manual(values=c("blue","orange"))+
	# theme_bw()+
	# labs(x="Lower rating",y="Higher rating",size="Count")+
	# theme(plot.margin = unit(c(5.5,-5.5,5.5,5.5), "pt"),axis.title.x=element_text(vjust=0,size=titles),axis.title.y=element_text(size=titles),
	# axis.text = element_text(size = txt),
	# legend.text = element_text(size = txt*.75),
	# legend.title = element_text(size=titles),
	# legend.box.margin=margin(0,0,0,-10))
