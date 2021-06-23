
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))
  # create parser object
parser <- ArgumentParser()

  # specify our desired options 
  # by default ArgumentParser will add an help option 
parser$add_argument("-s", "--snp_list", type="character", 
	help="tab sep 1-based SNP list file")
parser$add_argument("-tr", "--t_reference", type="character",
    help = "transcriptome reference fasta file")
parser$add_argument("-gr", "--g_reference", type="character",
    help = "genome reference fasta file")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

snp_list <- args$snp_list
tref <- args$t_reference
gref <- args$g_reference

#out_name <- gsub(".txt$", "", snp_list)

if(!file.exists(snp_list)){
  stop("no such snp list file exists")
}
if(!file.exists(tref)){
  stop("no such transcriptome reference file exists")
}
if(!file.exists(gref)){
  stop("no such genome reference file exists")
}

refseq <- readDNAStringSet(tref)
g_refseq <- readDNAStringSet(gref)

g_refseq_len <- tibble(g_CHROM=gsub(" .*", "", names(g_refseq)), g_len=width(g_refseq)) #genomic ref sequence lengths

snp_list <- read_table2(snp_list, col_names = TRUE)
colnames(snp_list)[1]="CHROM"
snp_list$SNP_Name <- paste0("S", str_pad(c(1:nrow(snp_list)), 8, pad="0"))
snp_list$CHROM_LEN <- width(refseq)[match(snp_list$CHROM, names(refseq))]


bed_left <- snp_list %>% mutate(start=POS-101, end=POS-1, name=paste0(SNP_Name, "-lf")) %>% select(CHROM, start, end, name, CHROM_LEN)
bed_left$start[bed_left$start<1]=1
bed_left <- bed_left %>% filter(end>=40) #keep only those where there is 40bp or more of left flank seq

bed_right <- snp_list %>% mutate(start=POS, end=POS+100, name=paste0(SNP_Name, "-rf")) %>% select(CHROM, start, end, name, CHROM_LEN)
bed_right$end[bed_right$end>bed_right$CHROM_LEN]=bed_right$CHROM_LEN[bed_right$end>bed_right$CHROM_LEN]
bed_right <- bed_right %>% filter(start<(CHROM_LEN-40)) #keep those where there is 40bp or more of right flank seq

bed <- bind_rows(bed_left, bed_right) %>% select(-CHROM_LEN) %>% arrange(name)
names(bed)[1]="chr"

write.table(bed, "flankseq.bed", row.names=F, col.names=F, quote=F, sep="\t")

snpContig_refseq <- refseq[names(refseq)%in%snp_list$CHROM]
writeXStringSet(snpContig_refseq, 'snpContig_refseq.fasta')

cat("\n Transcriptome SNPs contained within", length(snpContig_refseq), "contigs \n")

cat("Extracting flanking sequences\n")

command=paste("bedtools getfasta", "-fi", tref, "-bed", "flankseq.bed", "-name", "> flankseq.fasta", sep=" ")
try(system(command))

if(!file.exists("flankseq.fasta")){ stop("flanking seq fasta missing") }

cat("Running BLAST of flanking seq\n")

command=paste("makeblastdb", "-in", gref, "-dbtype", "nucl", "-parse_seqids", "-out", "refgenome.db.fasta", sep=" ")
try(system(command))

command=paste("blastn", "-query", "flankseq.fasta", "-db", "refgenome.db.fasta", "-outfmt", "6", "-evalue", "1e-6", "-out", "flankingSeq_blastn_output.txt", sep=" ")
try(system(command))

if(!file.exists("flankingSeq_blastn_output.txt")){ stop("flanking seq BLAST output missing") }

cat("Running BLAST of SNP contigs\n")

command=paste("blastn", "-query", "snpContig_refseq.fasta", "-db", "refgenome.db.fasta", "-outfmt", "6", "-evalue", "1e-6", "-out", "snpContig_blastn_output.txt", sep=" ")
try(system(command))

if(!file.exists("snpContig_blastn_output.txt")){ stop("SNP contig BLAST output missing") }


#blast reference transcript contig to genome to use to support flanking seq matches

blastn_flank <- read_table2("flankingSeq_blastn_output.txt", col_names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

blastn_snpContig <- read_table2("snpContig_blastn_output.txt", col_names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))



#take multi hits to same scaffold/contig and check how much of the snp contig is represented. 
t2g_coverage <- function(blastn_snpContig, snpContig_refseq, c, min.coverage=0.8){
	data <- blastn_snpContig %>% filter(qseqid==c)
	res <- NULL
	for(i in unique(data$sseqid)){
		qrange <- data %>% filter(sseqid==i) %>% select(qstart, qend)
		idx <- unique(unlist(apply(qrange, 1, function(x) seq(x[1], x[2]))))
		coverage <- prop.table(table(seq(1,width(snpContig_refseq[c]))%in%idx))["TRUE"]
		if(coverage>=min.coverage){
			res <- rbind(res, c(c, i, coverage))
		}
	}
	if(is.null(res)){
		res <- matrix(c(c, NA, NA),nrow=1)
	}
	res <- as_tibble(as.data.frame(res, stringsAsFactors=FALSE))
	names(res)=c("qseqid", "sseqid", "coverage")
	res$coverage <- as.numeric(res$coverage)	
	res
}

gregion <- bind_rows(lapply(unique(blastn_snpContig$qseqid), function(x) t2g_coverage(blastn_snpContig, snpContig_refseq, c=x, min.coverage=0.8)))

if(nrow(gregion)==0){ stop("No genome matches found for SNP contigs") }
cat("\n Transcriptome SNPs contained within", length(snpContig_refseq), "contigs found within", nrow(gregion), "genome scaffolds\n")



position_snp <- function(blastn_flank, gregion, snp_list, g_refseq, s, min.ident=0.85){
	data <- bind_rows(blastn_flank %>% filter(qseqid==paste0(s,"-lf"), qend==100, pident>=min.ident), blastn_flank %>% filter(qseqid==paste0(s,"-rf"), qstart==1, pident>=min.ident)) %>% mutate(g_region_match=sseqid%in%gregion$sseqid[gregion$qseqid==snp_list$CHROM[snp_list$SNP_Name==s]])
	res <- NULL
	if(nrow(data)!=0){
		for(i in unique(data$sseqid)){
			sdata <- data %>% filter(sseqid==i)
			if((sdata %>% count(qseqid) %>% pull(n) %>% max)==1){
				sdata <- sdata %>% mutate(pos=if_else(send > sstart, if_else(grepl("-lf", qseqid), send+1, sstart-1), if_else(grepl("-lf", qseqid), send-1, sstart+1))) # get snp pos factoring in that the genomic ref match may be in rev orientation to the transcript flanking seq
				if(nrow(sdata)==2){ # left and right flanked
						if(sdata$pos[1]==sdata$pos[2]){
							res <- rbind(res, c(s, i, sdata$pos[1]-1, sdata$pos[1], sdata$g_region_match[1],"dual_flanked"))
						}else{
							#error left and right flank do not match same pos
							res <- rbind(res, c(s, i, NA, NA, sdata$g_region_match[1], "left_right_mismatch"))
						}
				}else{ #single flanked
						if(grepl("-lf", sdata$qseqid[1])){ #left flanked
							res <- rbind(res, c(s, i, sdata$pos[1]-1, sdata$pos[1], sdata$g_region_match[1], "left_flanked"))
						}else{ #right flanked
							res <- rbind(res, c(s, i, sdata$pos[1]-1, sdata$pos[1], sdata$g_region_match[1], "right_flanked"))
						}
				}
			}else{
				#error multi hits within region
				res <- rbind(res, c(s, i, NA, NA, sdata$g_region_match[1], "multi_flanked_hits"))
			}
		}
		res <- as_tibble(as.data.frame(res, stringsAsFactors=FALSE))
		names(res)=c("SNP_Name", "g_CHROM", "bed_start", "bed_end", "snpContig_gregion_match", "Comment")
		res
	}
}


gSNP_pos <- bind_rows(lapply(unique(snp_list$SNP_Name), function(x) position_snp(blastn_flank, gregion, snp_list, s=x, min.ident=0.85)))

gSNP_pos_unfiltered <- gSNP_pos # store unfiltered gSNP_pos to allow manual filtering when loading in workspace image

gSNP_pos$bed_start <- as.numeric(gSNP_pos$bed_start)
gSNP_pos$bed_end <- as.numeric(gSNP_pos$bed_end)

gSNP_pos <- left_join(gSNP_pos, g_refseq_len, by="g_CHROM") %>% mutate(Comment = ifelse((bed_end>g_len | bed_end==0), "edge_overhang", Comment)) %>% select(-g_len)

if(nrow(gSNP_pos)==0){ stop("No genome positions found SNPs") }
cat("\n", gSNP_pos %>% n_distinct("SNP_Name"), "SNPs found with flanking genome hits out of a total", nrow(snp_list), "initial SNPs\n")


#filter for 2x locus for snp to report and remove non snpContig_gregion_match hits
gSNP_pos <- gSNP_pos %>% filter(snpContig_gregion_match==TRUE) %>% filter(Comment%in%c("dual_flanked", "left_flanked", "right_flanked")) %>% add_count(SNP_Name, name="g_hits") %>% filter(g_hits<=2)

if(nrow(gSNP_pos)==0){ stop("No genome SNPs positions after remove >2 g_hits") }
cat("\n", n_distinct(gSNP_pos$SNP_Name), "SNPs remaining after removing those with >2 genome hits or flanking errors\n")

gSNP_pos %>% distinct(SNP_Name, g_hits) %>% count(g_hits, name="No.SNPs")

gSNP_pos <- left_join(gSNP_pos, snp_list, by="SNP_Name") %>% select(-c(snpContig_gregion_match,CHROM_LEN))
colnames(gSNP_pos)[which(colnames(gSNP_pos)=="CHROM")] <- "t_CHROM"
colnames(gSNP_pos)[which(colnames(gSNP_pos)=="POS")] <- "t_POS"
colnames(gSNP_pos)[which(colnames(gSNP_pos)=="REF,ALT")] <- "t_REF,t_ALT"

write.table(gSNP_pos, "genomic_snpPos_info.txt", row.names=F, quote=F)

save.image("t2g-snp-workspace.Rdata")

