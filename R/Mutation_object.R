#' Add sample info to the MutationObject
#' @description This function is designed for merge multiple snp_anno and the
#' corresponding snp_data files, especially for multiple snp_anno and snp_data
#' generated from the Integrate_data() function.
#' @param message MutationObject.
#' @param verbose Sample info dataframe to add. Must have one column named
#' Sample_name which contain samples name identical to snp_data/snp_anno.
print_message <- function(message,verbose){
  if (verbose==T){
    message(message)
  }
}



#' @title The Mutation object for snp_data and snp_anno
#' @description The Mutation object is created to apply all functions for mutation analysis
#' @export
#'
SNP <- R6Class(
  classname = "Mutation Object",
  public = list(
    #' @field snp_list A numeric field
    #' @field metadata A numeric field
    #' @field modification_info A non-numeric field
    #' @field nj_distance NJ distance
    snp_list = NA,
    metadata = NULL,
    nj_distance = NA,
    modification_info = list(),
    #' Initialize a new MyClass object
    #' @param snp_anno Initial value for the field snp_anno
    #' @param snp_data Initial value for the field snp_data
    #' @param metadata Initial value for the field metadata
    #' @param remove_genes_dt Initial value for the field remove_genes_dt
    #' @param modification_info Initial value for the field remove_genes_dt
    #' @param annotate_AA Initial value for the field remove_genes_dt
    #' @param verbose Initial value for the field remove_genes_dt
    initialize = function(snp_anno,snp_data,metadata=NULL,
                          remove_genes_dt=NA,
                          annotate_AA=T,verbose=F) {

      self$metadata <- metadata

      ###### 1.sample info initialization

      if (!"REF" %in% colnames(snp_anno)){
        row_names <- rownames(snp_anno)
        row_names_unique_ind <- !duplicated(row_names)
        snp_anno <- snp_anno[row_names_unique_ind,]
        snp_data <- snp_data[,row_names_unique_ind,drop=F]
        rownames(snp_anno) <- row_names[row_names_unique_ind]
        colnames(snp_data) <- row_names[row_names_unique_ind]

      } else{
        row_names <- paste0(snp_anno$CHROM,"_",snp_anno$REF,snp_anno$POS,snp_anno$ALT,snp_anno$TYPE)
        rownames(snp_anno) <- row_names
        colnames(snp_data) <- row_names
        # snp_anno <- snp_anno %>% dplyr::select(-c("REF","ALT"))
      }

      snp_anno$INDEX <- rownames(snp_anno)

      # remove PE/PPE、phages, or repeated sequences
      if (is.na(remove_genes_dt)){
        print_message("[1/2] No remove_genes_dt, skip...",verbose)
        snp_anno_std <- snp_anno
      } else if (remove_genes_dt=="Mtb"){
        #data("Mtb_removed_gene_dt", package = "qMut")

        removed_genes <- Mtb_removed_gene_dt

        snp_anno_std <- snp_anno %>%
          dplyr::filter(!LOCUS_TAG %in% removed_genes$Locus)

        print_message("[1/2] Remove selected genes",verbose)

      } else {
        removed_genes <- fread(remove_genes_dt,data.table = F) %>%
          filter(Annotation!="")

        snp_anno_std <- snp_anno %>%
          dplyr::filter(!LOCUS_TAG %in% removed_genes$Locus)

        snp_anno_std <- snp_anno %>%
          dplyr::filter(!LOCUS_TAG %in% removed_genes$Locus)

        print_message("[1/2] Remove selected genes",verbose)

      }

      snp_anno_std <- snp_anno_std %>%
        dplyr::mutate(GENE=ifelse(is.na(GENE)|GENE=="",LOCUS_TAG,GENE))

      if (annotate_AA==T){
        print_message("[2/2] Annotate amino acid changes",verbose)

        #data("AA_3to1_dt", package = "qMut")


        AA_name <- AA_3to1_dt

        AA_effect <- str_extract(snp_anno_std$EFFECT,"(?<=p\\.).+")

        snp_anno_std$AA_effect_short <-
          stringi::stri_replace_all_fixed(AA_effect,AA_name$three,AA_name$one,vectorize_all = F)


      }

      snp_data_std <- snp_data[,rownames(snp_anno_std), drop=FALSE]

      self$snp_list <- list(snp_data=snp_data_std,snp_anno=snp_anno_std)

      # if (self$drug_name != "all"){
      #   selected_sample <- initial_phenotype$sample
      #
      #   if(sum(!initial_phenotype$sample %in% rownames(snp_data))>0){
      #     cat(paste0(sum(!initial_phenotype$sample %in% rownames(snp_data)),
      #                " sample were removed as they did not have mutation info!!!"))
      #
      #     selected_sample <- selected_sample[selected_sample %in% rownames(snp_data)]
      #   }
      #
      #   self$snp_list$snp_data <- self$snp_list$snp_data[initial_phenotype$sample,]
      #
      #
      #   self$snp_list$snp_data <- self$snp_list$snp_data[,Matrix::colSums(self$snp_list$snp_data)>0]
      #   self$snp_list$snp_anno <- self$snp_list$snp_anno[colnames(self$snp_list$snp_data),]
      #
      #   cat(paste("total",dim(self$snp_list$snp_data)[1],
      #             "sample and" ,dim(self$snp_list$snp_data)[2],"snps\n"))
      # } else {print_message("Choose all samples,so that all samples are remained",verbose)}

      if (class(self$metadata)=="data.frame"){
        rownames(metadata) <- metadata$Sample_name

        metadata <- data.frame(Sample_name=rownames(self$snp_list$snp_data)) %>%
          left_join(metadata,by="Sample_name")

        self$metadata <- metadata
      }

    },

    # update_obj = function(INDEX){
    #   new_obj <- self$clone(deep=T)
    #   new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[,INDEX,drop=F]
    #   new_obj$snp_list$snp_anno <- new_obj$snp_list$snp_anno[INDEX,]
    #   return(new_obj)
    # },

    #' get_metadata
    #' @param
    get_metadata = function(){

      if (!is.null(self$metadata) && (is.data.frame(self$metadata))){
        metadata_new <- self$metadata %>% filter(Sample_name %in% rownames(self$snp_list$snp_data))
        return(metadata_new)
      }

      return(NULL)


    },
    #' Initialize a new MyClass object
    #' @param selected_pheno Initial value for the field selected_pheno
    #' @return NULL
    ML_data_dt = function(selected_pheno){

      phenotype <- ifelse(self$phenotype$phenotype=="S",0,1)
      ML_data_dt <- cbind(phenotype,self$snp_list$snp_data)

      return(ML_data_dt)
    },

    #' Initialize a new MyClass object
    #' @param x Initial value for the field x
    Evolution_dt = function(x){
      #cat("[1/2]Melt data... \n")
      melt_dt <- mefa4::Melt(Matrix::t(self$snp_list$snp_data))

      #cat("[2/2]Create Evolution dataframe... \n")
      Evolution_dt <- data.frame(sampleID=melt_dt$cols,
                                 chr="NC_000962.3",
                                 pos=str_extract(melt_dt$rows,"\\d+"),
                                 ref=str_extract(melt_dt$rows,"[A-z]+"),
                                 mut=str_extract(melt_dt$rows,"(?<=\\d)[A-Z]+"))
      return(Evolution_dt)
      # mnt_line <- Evolution_dt %>% filter(nchar(ref)==nchar(mut)&nchar(ref)>1&nchar(mut)>1)
      # mnt2snp <- function(row_id){
      #   temp_line <- mnt_line[row_id ,]
      #   temp <- temp_line[c("ref","mut")]
      #   split_code <- str_split(temp,"")
      #   changed_code_index <- split_code[[1]]!=split_code[[2]]
      #   data.frame(
      #     sampleID = temp_line$sampleID,
      #     chr = temp_line$chr,
      #     pos = as.numeric(temp_line$pos)+c(1:length(changed_code_index)-1)[changed_code_index],
      #     ref = split_code[[1]][changed_code_index],
      #     mut = split_code[[2]][changed_code_index]
      #   )
      # }
      # mnt2snp_res <- pblapply(1:nrow(mnt_line),mnt2snp) %>% rbindlist() %>% as.data.frame()
      #
      # Evolution_dt <- Evolution_dt %>%
      #   filter(!(nchar(ref)==nchar(mut)&nchar(ref)>1&nchar(mut)>1)) %>%
      #   rbind(mnt2snp_res)



    },

    #' Getting SNPs in certain Genes, or certain SNP in all samples
    #' @param gene values for Index , Locus_tag or Gene
    #' @param gene_type values for searching mutations by certain type, should be one of
    #' "GENE","LOCUS_TAG","INDEX"
    #' @param upstream_length whether to show the upstream mutation of genes
    #' @param length Length of upstream and downstream regions of genes that should be included
    filter_snp = function(gene,gene_type="INDEX",upstream_length=NULL){ # must use =, or All elements of public, private, and active must be named.
      anno <- self$snp_list$snp_anno
      if (gene_type %in% c("GENE","LOCUS_TAG")){
        cal_promoter_region <- function(g,upstream_length){
          info <- anno %>%
            filter(LOCUS_TAG == !!g | GENE==!!g) %>%
            dplyr::slice(1)
          if (info$STRAND=="+"){
            gene_start <- info$POS-as.numeric(str_split(info$NT_POS,"/")[[1]][1])+1
            gene_stop <- info$POS-as.numeric(str_split(info$NT_POS,"/")[[1]][1])+as.numeric(str_split(info$NT_POS,"/")[[1]][2])

            promoter <- c(gene_start-upstream_length):c(gene_start-1)
          } else if (info$STRAND=="-"){
            gene_start <- info$POS+as.numeric(str_split(info$NT_POS,"/")[[1]][1])-1

            promoter <- c(gene_start+1):c(gene_start+upstream_length)
          } else {
            promoter <- c()
          }
          index_promotor <- self$snp_list$snp_anno %>%
            filter(POS %in% promoter) %>%
            rownames()
          if (length(index_promotor)==0){
            index_promotor <- NA
          }
          return(index_promotor)
        }
        index <- anno %>%
          filter(GENE %in% !!gene | LOCUS_TAG %in% !!gene) %>%
          rownames()

        if (!is.null(upstream_length)){
          index_promoter <- unlist(map2(gene,upstream_length,cal_promoter_region))
          index <- na.omit(unique(c(index,index_promoter)))
         }

      } else if (by=="INDEX"){index <- gene}

      else {
        stop("Parameter 'by' should be one of POS,GENE,LOCUS_TAG,INDEX!!!")
      }

      new_obj <- self$clone(deep=T)
      new_obj$snp_list$snp_anno <- new_obj$snp_list$snp_anno[index,]
      new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[,index,drop=F]

      new_obj$metadata <- new_obj$get_metadata()

      if (is.null(dim(new_obj$snp_list$snp_data))){
        cat(paste("remain",length(new_obj$snp_list$snp_data),
                  "sample and" ,"1","snps\n"))
      } else{
        cat(paste("remain",dim(new_obj$snp_list$snp_data)[1],
                  "sample and" ,dim(new_obj$snp_list$snp_data)[2],"snps\n"))
      }

      return(new_obj)

    }, # 筛选特定snp

    #' Initialize a new MyClass object
    #' @param sample_name sample names to filter
    #' @param mutation_index mutation index.
    #' @param opposite Used with mutation_index. When opposite is set to TRUE, samples
    #' containing mutation_index will be removed; otherwise will be remained.
    filter_sample = function(sample_name=NULL,mutation_index=NULL,opposite=F){

      new_obj <- self$clone(deep=T)

      if (!is.null(sample_name)){
        if (is.null(dim(new_obj$snp_list$snp_data))){
          new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[sample_name]

          new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[new_obj$snp_list$snp_data>0]
          new_obj$snp_list$snp_anno <- new_obj$snp_list$snp_anno[names(new_obj$snp_list$snp_data),]

          cat(paste("remain",length(new_obj$snp_list$snp_data),
                    "sample and" ,"1","snps\n"))
        }else{
          new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[sample_name,,drop=F]

          new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[,Matrix::colSums(new_obj$snp_list$snp_data)>0,drop=F]
          new_obj$snp_list$snp_anno <- new_obj$snp_list$snp_anno[colnames(new_obj$snp_list$snp_data),]

          cat(paste("remain",dim(new_obj$snp_list$snp_data)[1],
                    "sample and" ,dim(new_obj$snp_list$snp_data)[2],"snps\n"))
        }

      }

      if (!is.null(mutation_index)){
        # new_obj <- new_obj$filter_snp(mutation_index,"INDEX")

        info <- new_obj$snp_list$snp_data[,mutation_index]

        if (opposite==F){
          if (is.null(dim(info))){
            selected_samples <- names(info[info>0])
          } else{
            selected_samples <- rownames(info)[Matrix::rowSums(info)>0]
          }

        } else {
          if (is.null(dim(info))){
            selected_samples <- names(info[info==0])
          } else{
            selected_samples <- rownames(info)[Matrix::rowSums(info)==0]
          }
        }

        new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[selected_samples,,drop=F]
        # Matrix:colSums 和矩阵转置t()一样，都要用Matrix:，不然会因为调用base函数报错，因为无法处理稀疏矩阵
        new_obj$snp_list$snp_data <- new_obj$snp_list$snp_data[,Matrix::colSums(new_obj$snp_list$snp_data)>0]
        new_obj$snp_list$snp_anno <- new_obj$snp_list$snp_anno[colnames(new_obj$snp_list$snp_data),]

        cat(paste("remain",dim(new_obj$snp_list$snp_data)[1],
                  "sample and" ,dim(new_obj$snp_list$snp_data)[2],"snps\n"))
      }

      new_obj$metadata <- new_obj$get_metadata()

      return(new_obj)

    }, # 筛选特定样本，且仅保留至少存在一次的snp


    estimate_imp_all = function(){

      Times <- Matrix::colSums(self$snp_list$snp_data)
      return(
        self$snp_list$snp_anno %>%
          mutate(Times=Times)
      )
    },
    #' Getting SNPs counts in all samples
    #' @param gene values for Index , Locus_tag or Gene
    #' @param gene_type values for searching mutations by certain type, should be one of
    #' "GENE","LOCUS_TAG","INDEX"
    #' @param metadata_col Default is NULL. If the metadata_col was set, the estimate_imp function
    #' will calculate significance difference of mutation counts among different metadata_col levels.
    #' Should be colnames in metadata or the mutation INDEX.
    #' @param upstream_length whether to show the upstream mutation of genes
    estimate_imp = function(gene,gene_type="GENE", metadata_col=NULL,upstream_length=NULL){

      if (length(gene)==1 && gene=="all" && is.null(metadata_col)){

        estimate_res <- self$estimate_imp_all()
        return(estimate_res)

      }

      new_obj <- self$clone(deep=T)

      all_index <- colnames(self$snp_list$snp_data)

      if (gene!="all"){
        if (gene_type %in% c("GENE","LOCUS_TAG")){
          message("filter genes")

          new_obj <- new_obj$filter_snp(gene=gene,gene_type=gene_type,upstream_length=upstream_length)
        }
      }



      if (is.null(metadata_col)){
        return(new_obj$estimate_imp_all())
      }



      message("calculate mutation frequency")

      if (!is.null(metadata_col)){
        if (metadata_col %in% all_index){

          group <- self$snp_list$snp_data[rownames(new_obj$snp_list$snp_data),metadata_col]

          tmp_dt <- cbind(group,new_obj$snp_list$snp_data)

          total_count <- Matrix::colSums(tmp_dt[,-1])

          positive_count <- Matrix::colSums(tmp_dt[tmp_dt[,"group"]==1,-1])

          negative_count <- total_count-positive_count

          group_1 <- sum(group == 1)
          group_0 <- sum(group == 0)

          # 在R中的比例
          percent_1 <- positive_count/group_1 %>% round(5)
          percent_0 <- negative_count/group_0 %>% round(5)

          snp_summary <- new_obj$snp_list$snp_anno %>%
            mutate(count_total=total_count,
                   count_1=positive_count,
                   count_0=negative_count,
                   n1=group_1,
                   n0=group_0) %>%
            mutate(
              P_value = mapply(function(k1, k0, n1, n0) {
                tbl <- matrix(c(k1, n1 - k1, k0, n0 - k0), nrow = 2, byrow = TRUE)
                if (all(tbl >= 5)) {
                  chisq.test(tbl, correct = FALSE)$p.value
                } else {
                  fisher.test(tbl)$p.value
                }
              },
              k1 = count_1,
              k0 = count_0,
              n1 = n1,
              n0 = n0)
            ) %>%
            mutate(FDR = p.adjust(P_value, method = "fdr")) #%>%
            # mutate(P_value=ifelse(P_value < 0.001, sprintf("%.2e", P_value), round(P_value, 3)),
            #        FDR=ifelse(FDR < 0.001, sprintf("%.2e", FDR), round(FDR, 3)))


          return(snp_summary)

        } else if (metadata_col %in% colnames(new_obj$metadata)){


          # 原始稀疏矩阵（行 = 样本，列 = SNPs）
          snp_mat <- new_obj$snp_list$snp_data
          group <- new_obj$metadata[,metadata_col]
          # 获取所有分组名
          group_levels <- unique(group)

          # 初始化结果矩阵
          snp_names <- colnames(snp_mat)
          result_list <- lapply(group_levels, function(g) {
            idx <- which(group == g)
            counts <- Matrix::colSums(snp_mat[idx, , drop = FALSE])
            total <- length(idx)
            percents <- round(counts / total, 5)
            data.frame(percents)
          })

          # 合并所有组的结果
          mutation_summary <- do.call(cbind, result_list)
          colnames(mutation_summary) <- paste("percents",group_levels,sep="_")

          mutation_summary <- cbind(new_obj$snp_list$snp_anno,mutation_summary)

          return(mutation_summary)

        } else{
          stop(paste(metadata_col,"is neither exist in metadata nor mutation INDEX!"))
        }
      }

    }, # 评估snp位点重要性
    #' MutTree
    #'plot phylogenetic tree with heatmap of snps and lineage composition
    #'
    #' @param locus_tag the locus_tag of wanted gene
    #' @param limit the minimum frequence a sample need to have,defaults to 5
    #' @param size the number of top occured mutation types to plot,defaults to 20
    #' @param metadata_col Default is NULL. If the metadata_col was set, the estimate_imp function
    #' will calculate significance difference of mutation counts among different metadata_col levels.
    #' Should be colnames in metadata or the mutation INDEX.
    MutTree = function(locus_tag,limit=5,size=20,metadata_col=NULL){

      new_obj <- self$clone(deep=T)

      if (!is.null(metadata_col)){
        lineage <- new_obj$metadata[,c("Sample_name",metadata_col)]
      }

      data_clean <- function(locus_tag){

        gene_filtered <- new_obj$filter_snp(gene = locus_tag,gene_type = "LOCUS_TAG",upstream_length = 0)

        ref <- gene_filtered$estimate_imp("all") %>%
          filter(!grepl("^synonymous",.$EFFECT)) %>%
          dplyr::select(AA_effect_short) %>%
          rownames_to_column() %>%
          drop_na(AA_effect_short)

        rep_name <- filter(ref,duplicated(ref$AA_effect_short)) %>%.[2] %>%
          unlist() %>% unique()

        snp_all <- gene_filtered$snp_list$snp_data %>%
          .[,colnames(.) %in% ref$rowname] %>%
          `colnames<-`(ref$AA_effect_short)

        #对同种突变进行合并，再合并为所有重复突变的matrix与不含重复突变的矩阵合并
        snp_all <- map(rep_name,function(x){
          Matrix::rowSums(snp_all[,colnames(snp_all)==x]) %>%
            as.matrix() %>%
            `colnames<-`(x) %>% return()}) %>%
          do.call(cbind,.) %>%
          cbind(snp_all[,!(colnames(snp_all) %in% rep_name)])
        #将所有大于1的改为1
        snp_all[snp_all>1] <- 1

        cat(paste0("remain ",ncol(snp_all) ," significant snps\n"))

        return(snp_all)
      }

      snp_sort <- function(snp_all,limit=5){
        res_lineage <- NULL
        #单突变
        #通过@p相减得到每列突变样本数
        single <- snp_all@i[!duplicated(snp_all@i) & !duplicated(snp_all@i,fromLast=T)] + 1
        snp_single <- snp_all[single,]

        single_times <- diff(snp_single@p) %>%
          as.data.frame() %>% setNames("Freq") %>%
          rownames_to_column("j") %>%
          filter(Freq >= limit) %>%
          arrange(desc(Freq)) %>%
          dplyr::mutate(AA_effect_short=colnames(snp_single)[as.integer(.$j)],
                        Sample_name=map(as.integer(.$j),~match(1,snp_single[,.])) %>%
                          unlist() %>%
                          rownames(snp_single)[.]) %>%
          dplyr::select(AA_effect_short,Freq,Sample_name)

        if (!is.null(metadata_col)){
          rownames(snp_single) <- filter(lineage,Sample_name %in%
                                           rownames(snp_single))[,metadata_col]
        }
        #单突变谱系信息
        single_lineage <- map(single_times$AA_effect_short,function(x){
          rownames(snp_single)[snp_single[,x]==1] %>% table() %>%
            as.data.frame() %>%
            mutate(AA_effect_short=x)}) %>%
          do.call(rbind,.) %>%
          setNames(c(metadata_col,"Freq","AA_effect_short"))

        #多突变
        #转化为dgTMatrix根据@j按行找到每个样本发生的突变
        multiple <- unique(snp_all@i[duplicated(snp_all@i) | duplicated(snp_all@i,fromLast=T)] + 1)

        snp_multiple <- snp_all[multiple,]

        if(length(multiple)==1){
          #对Rv0100这种只有一个多突变的进行处理
          snp_multiple_chr <- which(snp_multiple==1) %>%
            names() %>% paste0(collapse = ",") %>%
            {data.frame(Sample_name=rownames(snp_all)[multiple],AA_effect_short=.)}

          colnames(snp_multiple_chr)[1] <- "Sample_name"
        }
        else{
          snp_multiple_chr <- apply(snp_multiple,1,function(x){
            which(x==1) %>% names() %>% paste0(collapse = ",")}) %>%
            data.frame() %>%
            rownames_to_column("Sample_name")

          colnames(snp_multiple_chr)[2] <- "AA_effect_short"

        #对Rv0287这种无多突变的进行处理
        if(length(multiple)!=0){
          multiple_times <- snp_multiple_chr$AA_effect_short %>%
            table() %>%
            as.data.frame() %>%
            `colnames<-`(c("AA_effect_short","Freq")) %>%
            left_join(snp_multiple_chr %>% filter(!duplicated(.$AA_effect_short)),
                      by="AA_effect_short") %>%
            filter(Freq >= limit) %>%
            arrange(desc(Freq))

          if (!is.null(metadata_col)){

            #多突变谱系
            snp_multiple_chr <- left_join(snp_multiple_chr,lineage,by="Sample_name")

            multiple_lineage <- snp_multiple_chr %>%
              group_by(AA_effect_short,!!sym(metadata_col)) %>%
              summarise(Freq=n(),.groups="drop")
          }




          #合并得到res_times和res_lineage
          res_times <- full_join(single_times,multiple_times,
                                 by = c("AA_effect_short", "Freq","Sample_name")) %>%
            arrange(desc(Freq))

          #无突变label及次数
          res_times <- snp_all[Matrix::rowSums(snp_all)==0,] %>% {
            data.frame(AA_effect_short="WT",
                       Freq=nrow(.),
                       Sample_name=rownames(.)[1])} %>%
            full_join(res_times,by = join_by(AA_effect_short, Freq,Sample_name))

          if (!is.null(metadata_col)){

            res_lineage <- full_join(multiple_lineage,single_lineage,
                                     by= c("AA_effect_short",metadata_col,"Freq"))

            res_lineage <- filter(lineage,Sample_name %in%
                                    rownames(snp_all[Matrix::rowSums(snp_all)==0,])) %>%
              dplyr::select(metadata_col) %>% table() %>% as.data.frame() %>%
              dplyr::mutate(AA_effect_short="WT") %>%
              full_join(res_lineage,.,by=c("AA_effect_short", "Freq",metadata_col))

          }

          cat(paste0("remain ",nrow(res_times)," samples with Freq >= ",limit,"\n"))
        }
        else{
          res_times <- single_times
        }

        return(list(res_times=res_times,res_group=res_lineage))
        }

      }

      plottree <- function(res_times,res_lineage,snp_all,size=20){

        plot_times <- slice_max(res_times[-1,],Freq,n=size-1,with_ties = F) %>%
          full_join(res_times[1,],by = c("AA_effect_short", "Freq","Sample_name"))

        plot_matrix <- snp_all[plot_times$Sample_name,]
        rownames(plot_matrix) <- (plot_times$AA_effect_short)
        plot_matrix <- plot_matrix[,Matrix::colSums(plot_matrix)!=0]


        tree <- plot_matrix %>%
          dist() %>%
          nj() %>% root("WT",edgelabel = T)

        tree$edge.length[tree$edge.length < 0.02] <- 0

        order <- fortify(tree) %>% filter(isTip=="TRUE") %>%
          dplyr::select(label,y) %>% arrange(y) %>%
          dplyr::select(AA_effect_short=label) %>%
          left_join(plot_times[,1:2],by="AA_effect_short")

        if (!is.null(res_lineage)){
          plot_lineage <- filter(res_lineage,AA_effect_short %in% order$AA_effect_short)

          lineage_count <- table(lineage[,metadata_col]) %>%
            data.frame() %>%
            setNames(c(metadata_col,"Count"))

          plot_lineage <- plot_lineage %>%
            full_join(lineage_count,by=metadata_col) %>%
            dplyr::mutate(Normalized_freq=Freq/Count,
                          AA_effect_short=factor(AA_effect_short,levels = order$AA_effect_short))
        }


        plot_heatmap <- reshape2::melt(plot_matrix %>% as.matrix())

        p1 <-ggtree(tree,size=0.8,col="darkgreen")+
          geom_tiplab(aes(label=""),color="#E74C3C",align = T)+
          geom_nodepoint(aes(subset=branch.length!=0),color="orange",alpha=1,size=0.8)+
          theme(plot.margin = margin(0,-10,0,0))

        p2 <- ggplot(reshape2::melt(plot_matrix %>% as.matrix()),
                     aes(x=Var2,y=Var1,fill = value))+
          geom_tile(color="white",width=1,height=1)+
          theme_minimal()+
          labs(x=NULL,y=NULL)+
          theme(axis.text.y = element_blank(),
                axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,
                                           size = 5,
                                           margin = margin(b=0)),
                legend.position = "none",
                plot.margin = margin(0,0,0,0))+
          scale_y_discrete(limit=order$AA_effect_short)+
          scale_fill_gradient(low = "#A6CEE3",high = "#4DAF4A")
        
        color_pool <- c("#1F78B4","#FF7F00","#E31A1C","#CAB2D6","#A6CEE3","#FB9A99",
                        "#33A02C","#cdc0b0","#B2DF8A","#0048ba","#7ac5cd",
                        "#8D4C6A","#377EB8","#419681","#4DAF4A","#727E76","#984EA3",
                        "#CB6651","#FFBF19","#FFFF33","#D2AA2D","#A65628","#CE6B73",
                        "#F781BF")

        if (!is.null(metadata_col)){
          p3 <-ggplot(data=plot_lineage,
                      aes(x=Normalized_freq,y=AA_effect_short,fill=!!sym(metadata_col)))+
            geom_bar(stat = "identity",position = "fill") +
            labs(x="Normalized counts",y=NULL)+
            theme_minimal()+
            scale_fill_manual(values=color_pool)+
            theme(axis.text.y = element_text(size = 8,face = "bold"),
                  axis.text.x = element_blank(),
                  panel.grid = element_blank(),
                  plot.margin = margin(t=0,r=-2,b=0,l=0),
                  legend.key.size = unit(0.23,"inch"),
                  legend.text = element_text(size = 6),
                  legend.title = element_text(size=9))+
            scale_y_discrete(position = "right",
                             labels = lapply(order$Freq, function(freq) bquote(italic("n") ~ "=" ~ .(freq))))
        }



        plot <- ggdraw()+
          draw_plot(p1, x = 0, y = 0, width = 0.22, height = 1, scale = 1)+
          draw_plot(p2, x = 0.22, y = 0, width = 0.22, height = 1, scale = 1)

        if (!is.null(metadata_col)){
          plot <- plot +
            draw_plot(p3, x = 0.42, y = 0, width = 0.5, height = 1, scale = 1)

        }

        return(plot)

      }

      snp_all <- data_clean(locus_tag)
      res <- snp_sort(snp_all=snp_all,limit)

      res_times <- res$res_times
      res_lineage <- res$res_group
      size <- ifelse(size > nrow(res_times),nrow(res_times),size)

      plot <- plottree(res_times=res_times,
                       res_lineage=res_lineage,
                       snp_all=snp_all,
                       size=size)

      # ggsave("tree.png",plot,dpi=600,
      #        width = 14.6*max(size,20)/20,
      #        height = 7.5*max(size,20)/20,
      #        unit="cm",bg = "white")

      print(plot)

      res_times <-res_times %>% dplyr::select(-Sample_name)

      return(list(plot=plot,
                  res_times=res_times,
                  res_lineage=res_lineage,
                  max_size=nrow(res_times)))
    },

    #' MutNumPlot
    #' plot mutation times and modification types in a single gene.
    #' @param locus_tag the locus_tag of wanted gene
    #' @param limit the minimum frequence a sample need to have,defaults to 1
    #' @param label_num the number of top frequent mutationsto label, defaults to 10
    #' @param metadata_col Default is NULL. If the metadata_col was set, the estimate_imp function
    #' @param show_synomonous Whether to show synomonous mutations.
    #' @param Modification Modification types to show.
    MutNumPlot = function(locus_tag,limit=1,label_num=10,
                          show_synomonous=T,
                          Modification=NULL){
      #self <- Mtb_mut
      # Modification <- c("Acetylation","Phosphorylation")
      plot_df <- function(locus_tag,limit,label_num,Modification){
        mut_dt <- self$estimate_imp("all")

        #不含修饰信息的数据框
        if (show_synomonous==F){

          mut_dt <- mut_dt %>%
            filter(!str_detect(EFFECT,"synonymous_variant")) %>%
            filter((LOCUS_TAG==locus_tag | GENE==locus_tag) & Times > limit) %>%
            mutate(Codon_Pos=as.numeric(str_extract(AA_POS,"[^/]+"))) %>%
            dplyr::select(Codon_Pos,AA_effect_short,Times) %>%
            rownames_to_column("INDEX")
        } else {
          mut_dt <- mut_dt %>%
            filter((LOCUS_TAG==locus_tag | GENE==locus_tag) & Times > limit) %>%
            mutate(Codon_Pos=as.numeric(str_extract(AA_POS,"[^/]+"))) %>%
            dplyr::select(Codon_Pos,AA_effect_short,Times) %>%
            rownames_to_column("INDEX")
        }

        #修饰信息
        if(!is.null(Modification) & sum(Modification %in% colnames(self$snp_list$snp_anno))>0){
          modification <- self$snp_list$snp_anno %>%
            filter(LOCUS_TAG==locus_tag | GENE==locus_tag) %>%
            mutate(AA=str_sub(AA_effect_short,1,1)) %>%
            dplyr::select(AA,!!Modification) %>%
            filter(!if_all(c(!!Modification), is.na)) %>%
            mutate(modification=Modification[apply(.,1,function(x){which(x==TRUE)-1})] %>%
                     paste0(.)) %>%
            rownames_to_column("INDEX") %>%
            dplyr::select(INDEX,modification)

          mut_dt <- left_join(mut_dt,modification,by="INDEX") %>%
            arrange(desc(Times))

        } else{
          mut_dt$modification=NA
          mut_dt <- mut_dt %>%
            arrange(desc(Times))
          }

        # pho_Y_svg <- system.file(".",c("pho_Y.svg"), package = "qMut")
        # pho_S_svg <- system.file(".",c("pho_S.svg"), package = "qMut")
        # pho_T_svg <- system.file(".",c("pho_T.svg"), package = "qMut")
        # ace_svg <- system.file(".",c("ace.svg"), package = "qMut")
        # mut_dt$modification <- mut_dt$modification %>%
        #   {ifelse(.=="S,Phosphorylation",pho_S_svg,
        #           ifelse(.=="T,Phosphorylation",pho_T_svg,
        #                  ifelse(.=="Y,Phosphorylation",pho_Y_svg,
        #                         ifelse(.=="K,Acetylation",ace_svg,NA))))}

        mut_dt[label_num+1:nrow(mut_dt),"AA_effect_short"]  <-
          mut_dt[label_num+1:nrow(mut_dt),"AA_effect_short"] %>%
          ifelse(!is.na(mut_dt[label_num+1:nrow(mut_dt),"modification"]),
                 .,NA)

        return(mut_dt)
      }

      mut_dt <- plot_df(locus_tag,limit,label_num,Modification)
      color_pool <- c("#1F78B4","#FF7F00","#E31A1C","#CAB2D6","#A6CEE3","#FB9A99",
                      "#33A02C","#cdc0b0","#B2DF8A","#0048ba","#7ac5cd",
                      "#8D4C6A","#377EB8","#419681","#4DAF4A","#727E76","#984EA3",
                      "#CB6651","#FFBF19","#FFFF33","#D2AA2D","#A65628","#CE6B73",
                      "#F781BF")
      ggplot(mut_dt,aes(x=Codon_Pos,y=Times))+
        geom_point(data=mut_dt %>% filter(!modification %in% Modification),
                   mapping=aes(x=Codon_Pos,y=Times,color=log10(Times)),
                   na.rm = T)+
        # geom_point(data=mut_dt %>% filter(!is.na(modification)),
        #            mapping=aes(x=Codon_Pos,y=Times),color="red",alpha=0.8,
        #            na.rm = T)+
        # geom_image(aes(image=modification),na.rm = T,size=0.025)+
        theme_bw()+
        guides(color="none")+
        scale_y_log10()+
        scale_color_gradientn(colors = turbo(n=100,alpha=.5))+
        geom_text_repel(data=mut_dt %>% filter(is.na(modification)),aes(label=AA_effect_short),na.rm = T)+
        new_scale_color() +
        geom_text_repel(aes(label=AA_effect_short,color=modification),na.rm = T)+
        scale_color_manual(values=color_pool,na.translate = FALSE)+
        theme(legend.position = "top",legend.box = "horizontal")+
        labs(x=paste0("Codon Position of ", locus_tag),color="Modification type")
    },

    #' calculate_NJ
    #'calculate sample distance for object based on Neighbor-Joining Algorithm.
    #'@param num_thread number of threads to use
    #'@param snp_range SNPs within this frequency range were included in the NJ (Neighbor-Joining) distance calculation. This range was used to accelerate computational speed.
    calculate_NJ = function(num_thread=1,snp_range = c(0.05,0.95)){

      nj_res <- cal_nj(obj=self,num.thread=num_thread,snp_range=snp_range)

      return(nj_res)
    }



  ),

  private = list(

    initial_sample = NA

  ),

  active = list(
    # update_metadata = function(){
    #
    #   if (!is.null(self$metadata)){
    #     metadata_new <- self$metadata %>% filter(Sample_name %in% rownames(self$snp_list$snp_data))
    #     return(metadata_new)
    #   }
    #
    #   return(NULL)
    #
#
#     }
  )
)

#' A SNP object
#' @description This function create a basic SNP object designed for downstream
#' SNPs.analysis.
#'
#' @param snp_anno Required. Generated from the Integrate_data() function.
#' @param snp_data Required. Generated from the Integrate_data() function.
#'
#' @param metadata Optional.Samples in metadata without corresponding mutation identification will be automatically removed.
#' @param remove_genes_dt Optional. A dataframe containing removed genes. The column
#' of the dataframe should be matched with the example(Mtb_removed_gene_dt). Currently,
#' if you are working for mutation analysis of Mycobacterium tuberculosis, you can directly use presetted
#' remove_genes_dt through providing value "Mtb".
#'
#' @export
#'
CreateSNPObj <- function(snp_anno,snp_data,annotate_AA=T,
                         metadata=NULL,remove_genes_dt=NA,verbose=F){


  obj <- SNP$new(snp_anno = snp_anno,
                 snp_data = snp_data,
                 metadata = metadata,
                 annotate_AA = annotate_AA,
                 remove_genes_dt = remove_genes_dt,
                 verbose=verbose)
  return(obj)
}


create_counter_func <- function(obj_count) {
  count <- 0
  function() {
    count <<- count + 1
    cat("[",count,"/", obj_count,"]", "\n")

  }
}

merge_data_single <- function(SNP_obj1,SNP_obj2){

  # Build the integrated SNP_anno file
  cat("Merge snp_anno file \n")
  all_anno_dup <- rbind(SNP_obj1$snp_list$snp_anno %>% rownames_to_column("index"),
                        SNP_obj2$snp_list$snp_anno %>% rownames_to_column("index"))

  snp_anno_new <-  all_anno_dup %>%
    distinct(index,.keep_all = T) %>%
    arrange(POS) %>%
    column_to_rownames("index")

  # Build the integrated snp_data file
  cat("Merge snp_data file \n")
  anno_ind1 <- rownames(SNP_obj1$snp_list$snp_anno)
  anno_ind2 <- rownames(SNP_obj2$snp_list$snp_anno)

  anno1_anno2_shared <- intersect(anno_ind1,anno_ind2)
  anno1_unique <- setdiff(anno_ind1,anno_ind2)
  anno2_unique <- setdiff(anno_ind2,anno_ind1)

  data1 <- SNP_obj1$snp_list$snp_data
  data2 <- SNP_obj2$snp_list$snp_data

  # Extract common SNP data and integrate it
  data1_shared <- data1[,anno1_anno2_shared,drop=F]
  data2_shared <- data2[,anno1_anno2_shared,drop=F]

  data_shared <- rbind(data1_shared,data2_shared)

  # Retrieve the unique SNP data for Data 1, and at this point, all SNPs for Data 2 are marked as 0
  data1_unique <- data1[,anno1_unique,drop=F]
  data2_not_contain <- Matrix(0,
                              nrow=length(rownames(data2)),
                              ncol=length(anno1_unique),
                              dimnames = list(rownames(data2),anno1_unique))

  data_add1 <- rbind(data1_unique,data2_not_contain)

  # Retrieve the unique SNP data for Data 2, and at this point, all SNPs for Data 1 are marked as 0
  data2_unique <- data2[,anno2_unique,drop=F]
  data1_not_contain <- Matrix(0,
                              nrow=length(rownames(data1)),
                              ncol=length(anno2_unique),
                              dimnames = list(rownames(data1),anno2_unique))

  data_add2 <- Matrix::rbind2(data2_unique,data1_not_contain)

  # Merge the final new snp_data
  snp_data_new <- cbind(data_shared,data_add1,data_add2)[,rownames(snp_anno_new),drop=F]

  # Build the integrated metadata file
  cat("Merge metadata file \n")

  metadata1 <- SNP_obj1$metadata
  metadata2 <- SNP_obj2$metadata

  # If any sample information is missing, simply do not include the sample information
  if (is.null(metadata1) && is.null(metadata2)){
    metadata <- NULL
  } else {
    metadata <- plyr::rbind.fill(metadata1,metadata2)
  }

  all_data_merged <- SNP$new(snp_anno = snp_anno_new,
                             snp_data = snp_data_new,
                             metadata = metadata,
                             annotate_AA = F)

  # all_data_merged <- CreateSNPObj(snp_anno = snp_anno_new,
  #                                 snp_data = snp_data_new,
  #                                 metadata = metadata)

  return(all_data_merged)

}

#' Merge multiple snp_anno and snp_data files
#' @description This function is designed for merge multiple snp_anno and the
#' SNP Object.
#' @param obj_list Object list.
#' @export
MergeMutationObj <- function(obj_list){

  message("Total iteration time: ", length(obj_list)-1)
  batch_merge <- function(obj_list, batch_size = 2) {
    # 将数据框列表分成每 batch_size 一组
    batches <- split(obj_list, ceiling(seq_along(obj_list) / batch_size))

    # 对每一组进行合并
    merged_batches <- map(batches, function(batch) {
      purrr::reduce(batch, ~merge_data_single(.x,.y))
    },.progress = T)
    return(merged_batches)
}

  while (length(obj_list)>=2){

    message("Number of files is ", length(obj_list), ", combining files to accelerate the procedure")

    obj_list <- batch_merge(obj_list)
  }

  return(obj_list)
}


#' Add metadata to the MutationObject
#' @description This function is designed for add additional sample info for the object.
#' @param obj MutationObject.
#' @param metadata Sample info dataframe to add. Must have one column named
#' Sample_name which contain samples name identical to snp_data/snp_anno.
#' @export
AddMetadata <- function(obj,metadata){
  if (is.na(obj$metadata)){
    metadata_raw <- data.frame(Sample_name = rownames(obj$snp_list$snp_data))
  }
  modification_info <- left_join(metadata_raw,metadata,by="Sample_name")

  obj$metadata <- modification_info

  return(obj)

}



#' Add sample info to the MutationObject
#' @description This function is designed for merge multiple snp_anno and the
#' corresponding snp_data files, especially for multiple snp_anno and snp_data
#' generated from the Integrate_data() function.
#' @param obj SNPObject.
#' @param modification_dt Sample info dataframe to add. Must have one column named
#' @param type Sample info dataframe to add. Must have one column named
#' Sample_name which contain samples name identical to snp_data/snp_anno.
#' @export
AddModificationInfo <- function(obj,modification_dt,
                                type){

  message("Add ",type," info to MutationObject")
  # obj=Mtb_mut modification_dt=ace_info
  snp_anno_new <- obj$snp_list$snp_anno %>%
    mutate(AA_POS2=as.numeric(str_extract(AA_POS,"[^/]+")),
           INDEX=rownames(.))

  modification_dt[,type] <- T

  if (!"AA_POS" %in% colnames(modification_dt)){
    stop("Must provide the AA_POS column and GENE/LOCUS_TAG column!")
  }

  modification_dt <- modification_dt %>%
    dplyr::select(any_of(c("INDEX","AA_POS", "LOCUS_TAG","GENE","NT_POS",!!type))) %>%
    dplyr::rename(AA_POS2=AA_POS)

  # modification at amino acid level

  snp_anno_new <- suppressWarnings(
    left_join(snp_anno_new,modification_dt,
              by=intersect(c("LOCUS_TAG","GENE","NT_POS","AA_POS2"),colnames(modification_dt)))
  )
  snp_anno_new <- snp_anno_new %>%
    dplyr::select(-AA_POS2)


  # remove synonymous_variant
  # modification_info_new[[type]][str_detect(modification_info_new$EFFECT,"synonymous_variant")] <- NA

  # Sometimes, due to differences in the starting sites, some annotations may be incorrect. To be conservative, all were removed.
  if (type=="Phosphorylation"){
    snp_anno_new[[type]][snp_anno_new[[type]]==T &
                           !grepl("T|S|Y",str_extract(snp_anno_new$AA_effect_short,"[^\\d]+"))] <- NA
  } else if (type=="Acetylation"){
    snp_anno_new[[type]][snp_anno_new[[type]]==T &
                           !grepl("K",str_extract(snp_anno_new$AA_effect_short,"[^\\d]+"))] <- NA
  }

  snp_anno_new <- snp_anno_new %>%
    distinct(INDEX,.keep_all = T)

  rownames(snp_anno_new) <- snp_anno_new$INDEX
  obj$snp_list$snp_anno <- snp_anno_new

  return(obj)

}


