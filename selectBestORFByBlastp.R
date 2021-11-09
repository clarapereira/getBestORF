library(tidyverse)

basedir <- "/Users/clarapereira/Dropbox/Barreto_lab/ORFs/"

orfs.file <- paste0(basedir, "/orfs_out/orfs.tab")

peptides.path  <- paste0(basedir, "/orfs_out/orfipy_gencode.vM27.transcripts.fasta_out/")
peptides.file  <- paste0(peptides.path, "/orfs_peptides.tab")

blastp.path <- paste0(basedir, "/orfs_out/orfipy_gencode.vM27.transcripts.fasta_out/fasta/")
out.path <- paste0(basedir, "/results/")
                 
                      
                      
# 1. read in the ORFs

orfs.x <- data.table::fread(orfs.file)
orfs.pc.x <- orfs.x %>% #head()
  filter(V8 == "protein_coding"); orfs.pc.x %>% head()
rm(orfs.x)

orfs.pc <- orfs.pc.x %>% #head() %>% 
  separate(V9, into = c("ORFid", "ORFn"), sep = "[\t]") %>% 
  unite(col = "fullORFid", c(V1:ORFid), sep = "|", remove = F) %>% 
  unite(col = "shortORFid", c(V5,ORFid), sep = "", remove = F) %>% 
  select(-V1,-V2,-V3,-V4,-V5,-V8); orfs.pc %>% head()

rm(orfs.pc.x)

orfs.pc %>% nrow()


genes_list.pc <- orfs.pc %>% select(V6) %>% distinct() %>% unlist() %>% as.vector() 
length(genes_list.pc)
genes_list.pc %>% head(10)

#genes_list <- c("Gm21287", "Adhfe1")
#genes_list <- c("Gm21287", "Gm21287")

# data in chunks: 
chunk <- function(x,n){
  split(x, factor(sort(rank(x)%%n)))
} 

genes_list.pc.chunks <- chunk(genes_list.pc,1000)


# for (c in seq_along(genes_list.pc.chunks)){
#   genes_list <- genes_list.pc.chunks[[c]] 
#   print(genes_list)
#   path <-  paste0(out.path, "genes_list.", c,".tsv")
#   print(path)
#   write(genes_list, file = path)
# }

# Run the blast and ORF select operations over the chunks:
for (c in seq_along(genes_list.pc.chunks)){
  message <-  paste0("Staring chunk nr:  ", c)
  print(message)
  
  genes_list <- genes_list.pc.chunks[[c]]
  
  orfs_gene_best <- list()
  for (i in seq_along(genes_list)){
    
    # 0.1. extract the gene from the peptide tab file and convert it to .fasta
    # 0.2. run blastp
    # 0.3. do the R thing: 1, 2, 3
    gene_tab <- paste0(peptides.path, "/orfs_peptides.simpleIDs_",genes_list[i],".tab")
    gene_fasta <- paste0(peptides.path, "/orfs_peptides.simpleIDs_",genes_list[i],".fasta")
    extract_gene_from_tab <- paste0("grep ", genes_list[i]," ", peptides.file, " > ", gene_tab)
    system(
      # 3.4. Separate query .fasta file into N .fasta files (as many as the sequences)
      extract_gene_from_tab, 
      timeout = 300
    )
    
    tabtofasta <- paste0("python3 /Users/clarapereira/Dropbox/Barreto_lab/ORFs/convertTabToFasta.py ","", gene_tab, " ", gene_fasta)
    system(
      tabtofasta,
      #"pip3 install biopython",
      timeout = 500
    )
    
    splitfile <- paste0("sh fastaToManyFasta.sh ","", gene_fasta, " ", blastp.path)
    system(
      splitfile,
      timeout = 500
    )
    
    runblastp <- paste0("sh blastpOverFiles.sh ","", blastp.path, " ", "/Users/clarapereira/Dropbox/Barreto_lab/ORFs/reference/UP000000589_10090.fasta")
    
    system(
      runblastp,
      timeout = 1000
    )
    
    # 1. extract ORFs for a single gene 
    gene <- genes_list[i]
    orfs_gene <- orfs.pc %>% filter(V6 == gene) 
    print(gene)
    
    # 2. read in the blastp output and filter for the best match per ORF:
    blastp_gene.list <- list()
    for (n in seq_along(orfs_gene$fullORFid)){
      file <- paste0(blastp.path, orfs_gene$fullORFid[n], ".fa_blastpResult.fasta.tsv")
      
      blastp_gene.list[[n]] <- tryCatch({
        if (file.size(file) > 0){
          data.table::fread(file) %>% 
            filter( V11 == min(V11) & V12 == max(V12))
        }
      }, error = function(err) {
        # error handler picks up where error was generated
        print(paste("Read.table didn't work!:  ",err))
      }
      )
      
      #blastp_gene.list[[i]] <- data.table::fread(paste0(blastp.path, orfs_gene$fullORFid[i], ".fa_blastpResult.fasta.tsv")) %>% 
      # filter( V11 == min(V11) & V12 == max(V12))
      
    }
    # 2.1 merge all the best hits per ORF
    blastp_gene.list <- blastp_gene.list[blastp_gene.list != "Read.table didn't work!:   Error in if (file.size(file) > 0) {: missing value where TRUE/FALSE needed\n"]
    blastp_gene <- reduce(blastp_gene.list, rbind, fill=TRUE) 
    
    # output format description: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    # might be interesting also: https://www.biostars.org/p/88944/
    colnames(blastp_gene) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    
    # 3. Select the best ORF per gene: 
    
    orfs_gene_best[[i]] <- orfs_gene %>% 
      left_join(
        blastp_gene, 
        by = c("fullORFid" = "qseqid")
      ) %>% 
      na.omit() %>% 
      filter( evalue == min(evalue) & bitscore == max(bitscore) ) %>% #& pident == max(pident)
      distinct()
    
    message(
      "\n Number of best matches: ", orfs_gene_best[[i]] %>% nrow(),
      "\n********************************"
    )
    
    # 4. commandline: remove the data from the folder
    
    
    removefiles1 <- paste0("rm ", blastp.path, "*.fa")
    removefiles2 <- paste0("rm ", blastp.path, "*.fasta.tsv")
    removefiles3 <- paste0("rm ", gene_tab)
    removefiles4 <- paste0("rm ", gene_fasta)
    system(
      removefiles1,
      timeout = 500
    )
    system(
      removefiles2,
      timeout = 500
    )
    system(
      removefiles3,
      timeout = 500
    )
    system(
      removefiles4,
      timeout = 500
    )
  }
  length(orfs_gene_best)
  orfs_gene_best_merged <- reduce(orfs_gene_best, rbind, fill=TRUE) 
  
  path <- paste0(out.path, "bestORF_protein_coding.", c,".tsv")
  print(path)
  orfs_gene_best_merged %>% 
    write_delim(
      path,
      delim = "\t"
    )
}
# 
# 
# 
# 
# orfs_gene_best <- list()
# for (i in seq_along(genes_list)){
#   
#   # 0.1. extract the gene from the peptide tab file and convert it to .fasta
#   # 0.2. run blastp
#   # 0.3. do the R thing: 1, 2, 3
#   gene_tab <- paste0(peptides.path, "/orfs_peptides.simpleIDs_",genes_list[i],".tab")
#   gene_fasta <- paste0(peptides.path, "/orfs_peptides.simpleIDs_",genes_list[i],".fasta")
#   extract_gene_from_tab <- paste0("grep ", genes_list[i]," ", peptides.file, " > ", gene_tab)
#   system(
#     # 3.4. Separate query .fasta file into N .fasta files (as many as the sequences)
#     extract_gene_from_tab, 
#     timeout = 300
#     )
#   
#   tabtofasta <- paste0("python3 convertTabToFasta.py ","", gene_tab, " ", gene_fasta)
#   system(
#     tabtofasta,
#     timeout = 500
#   )
#   
#   splitfile <- paste0("sh fastaToManyFasta.sh ","", gene_fasta, " ", blastp.path)
#   system(
#     splitfile,
#     timeout = 500
#   )
#   
#   runblastp <- paste0("sh blastpOverFiles.sh ","", blastp.path, " ", "/Users/clarapereira/Dropbox/Barreto_lab/ORFs/reference/UP000000589_10090.fasta")
# 
#   system(
#     runblastp,
#     timeout = 1000
#   )
#   
#   # 1. extract ORFs for a single gene 
#   gene <- genes_list[i]
#   orfs_gene <- orfs %>% filter(V6 == gene) 
#   print(gene)
#   
#   # 2. read in the blastp output and filter for the best match per ORF:
#   blastp_gene.list <- list()
#   for (n in seq_along(orfs_gene$fullORFid)){
#     file <- paste0(blastp.path, orfs_gene$fullORFid[n], ".fa_blastpResult.fasta.tsv")
#     
#     blastp_gene.list[[n]] <- tryCatch({
#       if (file.size(file) > 0){
#         data.table::fread(file) %>% 
#           filter( V11 == min(V11) & V12 == max(V12))
#       }
#     }, error = function(err) {
#       # error handler picks up where error was generated
#       print(paste("Read.table didn't work!:  ",err))
#     }
#     )
#     
#     #blastp_gene.list[[i]] <- data.table::fread(paste0(blastp.path, orfs_gene$fullORFid[i], ".fa_blastpResult.fasta.tsv")) %>% 
#     # filter( V11 == min(V11) & V12 == max(V12))
#   
#   }
#   # 2.1 merge all the best hits per ORF
#   blastp_gene.list <- blastp_gene.list[blastp_gene.list != "Read.table didn't work!:   Error in if (file.size(file) > 0) {: missing value where TRUE/FALSE needed\n"]
#   blastp_gene <- reduce(blastp_gene.list, rbind, fill=TRUE) 
#   
#   # output format description: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
#   # might be interesting also: https://www.biostars.org/p/88944/
#   colnames(blastp_gene) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
#                              "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
#   
#   # 3. Select the best ORF per gene: 
#   
#   orfs_gene_best[[i]] <- orfs_gene %>% 
#     left_join(
#       blastp_gene, 
#       by = c("fullORFid" = "qseqid")
#     ) %>% 
#     na.omit() %>% 
#     filter( evalue == min(evalue) & bitscore == max(bitscore) ) %>% #& pident == max(pident)
#     distinct()
#   
#   message(
#     "\n Number of best matches: ", orfs_gene_best[[i]] %>% nrow(),
#     "\n********************************"
#   )
#   
#   # 4. commandline: remove the data from the folder
#   
#   
#   removefiles1 <- paste0("rm ", blastp.path, "*.fa")
#   removefiles2 <- paste0("rm ", blastp.path, "*.fasta.tsv")
#   system(
#     removefiles1,
#     timeout = 500
#   )
#   system(
#     removefiles2,
#     timeout = 500
#   )
# }
# length(orfs_gene_best)
# orfs_gene_best_merged <- reduce(orfs_gene_best, rbind, fill=TRUE) 
#   
# 
# orfs_gene_best_merged %>% 
#   write_delim(
#     paste0(out.path, "bestORF_protein_coding.", "1-10",".tsv"),
#     delim = "\t"
#   )
# 
# 
# 
# 
# 
# 
