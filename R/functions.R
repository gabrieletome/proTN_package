################################################################################
# ProTN: an integrative pipeline for complete analysis of proteomics           # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# Issue at: https://github.com/TebaldiLab/ProTN/issues                         #
# PI: Dr. Toma Tebaldi, PhD                                                    #
################################################################################

### Shared variables to load ----
colour_vec<-c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A","#FF7F00",
              "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F",
              "#9af7c9","#c0a0e6","#f5939a","#5fefee","#e8d388","#8abdff","#e3bf7a","#48c0f6","#ccb86f","#fafeaf",
              "#04bfe4","#ea92a8","#25cfe9","#ffbb97","#62e8ff","#ffd29b","#8fcbff","#aec87f","#b2c6ff","#aedf9b",
              "#ffb3c8","#81b97c","#ffcef0","#54bcae","#ffcab3","#7dd7ff","#ffe4ad","#81b0da","#ebffc5","#a7a8d5",
              "#bdffd0","#d6d2ff","#9fb282","#a5d8ff","#c9a395","#9bfff7","#ffcbce","#70b6ca","#faffdd","#88b3bc",
              "#d6ffdf","#a3e5ff","#b1ac8f","#ccffff","#b4bb9f","#7eb7a5","#c5dbc3","#a8b9ff","#74dfbc","#f2ace7")


### The function enrichment_test performs enrichment analysis starting from 3 vectors of gene names: input, annotation and background (the background is intended as the universe and, by default, is also used to filter gene names in the input and annotation) ----  


enrichment_test <- function(input_vec, # character vector, input genes
                            anno_vec, # character vector, annotation genes
                            background_vec, # character vector, background/universe genes
                            input_name="i1", # name for input list (i1 by default)
                            anno_name="a1", # name for annotation list (a1 by default)
                            anno_class="a", # name for annotation class (groups of annotations, a by default)
                            id_separator=";" # character separating overlap ids in the resulting string (; by default)
){
  bg_vec <- background_vec %>% str_trim() %>% str_to_upper() %>% 
    unique() %>% str_sort() # clean background
  
  ib_vec <- input_vec %>% str_trim() %>% str_to_upper() %>% 
    unique() %>% str_sort() %>% 
    extract(. %in% bg_vec) # clean input (overlap with bg)
  
  ab_vec <- anno_vec %>% str_trim() %>% str_to_upper() %>% 
    unique() %>% str_sort() %>% 
    extract(. %in% bg_vec) # clean anno (overlap with bg)
  
  ov_vec <- ib_vec %>% extract(. %in% ab_vec) # overlap between input and anno
  
  ovinput_vec <- input_vec %>% str_trim() %>% str_sort() %>% 
    unique() %>% extract(str_to_upper(.) %in% ov_vec) # overlap, with input case
  ovinput_vec <- ovinput_vec[!duplicated(str_to_upper(ovinput_vec))] # remove of case sensitive duplicates
  
  IY<- length(ov_vec) # in input, in anno (overlap)
  IN<- length(ib_vec) - IY # in input, not in anno
  BY<- length(ab_vec) - IY # not in input, in anno
  BN<- length(bg_vec) - (IY+IN+BY) # not in input, not in anno   
  
  # fisher test, 1 tail (only enrichment)
  fisher_data <- fisher.test(matrix(c(IY, BY, IN, BN), 2, 2), alternative="greater") 
  # hypergeometric test (equivalent to Fisher 1 tail)
  #phyper(q=(IY-1),m=(BY+IY),n=(IN+BN),k=(IY+IN),lower.tail = F) 
  
  # generation of the output (a tibble with 1 row)
  out_df <- tibble("input_name"=input_name, # name of the input
                   "anno_name"=anno_name, # name of the annotation
                   "anno_class"=anno_class,
                   "overlap_size"=IY, # size of the overlap 
                   "p_value"=fisher_data$p.value, # p-value of the enrichment (Fisher test)
                   "odds_ratio"=fisher_data$estimate, # conditional MLE estimate from fisher.test function
                   "combined_score"=0, # default combined score (to match enrichr enrichment results)
                   #"odds_ratio_sample"= (IY/IN)/(BY/BN), # sample odds ratio
                   #"fold_enrichment"= (IY/(IY+IN))/((IY+BY)/(BY+IY+BY+BN)), # fold enrichment formula
                   "input_size"=IY+IN, # size of processed input 
                   "anno_size"=IY+BY, # size of processed annotation
                   "background_size"= IY+IN+BY+BN, # size of processed background
                   "overlap_input_ratio"= IY/(IY+IN), # overlap ratio with respect to input size
                   "overlap_anno_ratio"= IY/(IY+BY), # overlap ratio with respect to annotation size
                   "overlap_ids"= str_c(ovinput_vec,collapse=id_separator) # string with overlap ids
  )
  return(out_df)
}



### The function enrichment_enrichr performs enrichment analysis with the enrichR bioconductor library, interface to the enrichR webserver ---- 

enrichment_enrichr <- function(input_vec, # vector with gene names
                               input_name="i1", # name for input list
                               dbs_vec=NULL, # vector with database names
                               back_frg = NULL,
                               return_size_in_anno = T # use only genes with at least 1 annotation to calculate input size
){ 
  
  input_vec <- input_vec %>% str_trim() %>% str_sort() %>% 
    unique() # clean input vector
  input_vec <- input_vec[!duplicated(str_to_upper(input_vec))] # remove case duplicates
  
  if(!is.null(back_frg)){
    back_frg <- back_frg %>% str_trim() %>% str_sort() %>% 
      unique() # clean input vector
    back_frg <- back_frg[!duplicated(str_to_upper(back_frg))] # remove case duplicates
  }
  
  dbs_default <- as.vector(suppressMessages(read_tsv(system.file("extdata", "dbs_enrichR.txt", package = "proTN"), col_names = FALSE))[,1]) %>% unlist()
  dbs_used <- tryCatch({
    dbs <- listEnrichrDbs()
    if(!is.null(dbs_vec)){dbs_used<-intersect(dbs_vec,dbs$libraryName)
    } else {dbs_used<-intersect(dbs_default,dbs$libraryName)}
    dbs_used
  }, error = function(cond){
    if(!is.null(dbs_vec)){dbs_used<-dbs_vec
    } else {dbs_used<-dbs_default}
    dbs_used
  })
  
  
  if(is.null(back_frg)){
    enrich_list <- enrichr(input_vec, dbs_used)
  } else {
    enrich_list <- enrichr(input_vec, dbs_used, background = back_frg, include_overlap = TRUE)
  }
  out_df <- ldply(enrich_list,.id="anno_class")
  rm(enrich_list)
  
  if(nrow(out_df) > 0){
  out_df <- out_df %>% rename(anno_name=Term,
                              p_value=P.value,
                              fdr=Adjusted.P.value,
                              odds_ratio=Odds.Ratio,
                              combined_score=Combined.Score,
                              enrichr_ids=Genes)  # rename columns
  
  transient <- str_split(out_df$Overlap,"/",simplify=T) # extract overlap and anno size
  out_df$overlap_size <- transient[,1] %>% as.numeric()
  out_df$anno_size <- transient[,2] %>% as.numeric()
  rm(transient)
  
  out_df <- out_df %>% select(-c(Old.P.value,Old.Adjusted.P.value,Overlap)) # remove useless columns
  out_df$overlap_anno_ratio <- out_df$overlap_size/out_df$anno_size
  
  # sort genes and recover original gene names used for input
  transient <- out_df %>% select(anno_class,anno_name,enrichr_ids) %>% 
    separate_rows(enrichr_ids,sep=";") %>% arrange(anno_class,anno_name,enrichr_ids)
  
  out_df$enrichr_ids<-NULL 
  out_df$input_size <- ifelse(return_size_in_anno, 
                              length(unique(transient$enrichr_ids)), 
                              length(input_vec)) # calculate input size
  
  match_df <- data.frame("input_ids"=input_vec,
                         "enrichr_ids"=toupper(input_vec))
  transient<- left_join(transient,match_df)
  
  transient <- transient %>% mutate(overlap_ids = ifelse(is.na(input_ids),enrichr_ids,input_ids))
  transient <- transient %>% group_by(anno_class,anno_name) %>% 
    summarise("overlap_ids"=str_c(overlap_ids,collapse=";")) %>% ungroup() 
  
  out_df <- left_join(out_df,transient)
  rm(transient,match_df)
  
  out_df$overlap_input_ratio<-out_df$overlap_size/out_df$input_size
  out_df$input_name<- input_name
  
  out_df <-out_df %>% select(input_name,anno_name,anno_class,
                                    overlap_size,p_value,fdr,odds_ratio,combined_score,
                                    input_size,anno_size,overlap_input_ratio,overlap_anno_ratio,
                                    overlap_ids) %>% as_tibble()
  # background_size is not estimated, odds ratio retrived from enrichr
  } else{
    out_df <- NULL
  }
  
  return(out_df)
}



### The function enrichment_lollipop creates an enrichment plot (lollipop) starting from an enrichment dataframe ---- 

enrichment_lollipop <- function(input_df, # input dataframe
                                x_col="combined_score", # column for x axis, necessary
                                y_col="anno_name", # column for y axis, necessary
                                
                                size_col="overlap_size", # column for size, optional
                                size_vec=c(1,5), # size vector c(min,max), optional
                                
                                shape_col="P<0.05", # column for shape, optional
                                shape_vec=c(16,21), # shape vector, optional
                                
                                color_col="P<0.05", # column for color, optional
                                color_vec=c("grey25","grey50"), # color vector, optional
                                
                                fill_col="P<0.05", # column for fill color, optional
                                fill_vec=c("grey25","white"), # fill vector, optional
                                
                                text_col="overlap_size", # column for fill color, optional
                                
                                facet_col="anno_class", # column for facet, optional
                                size_text=3,
                                char_max = 50, # maximum number of chars in text
                                sort_y = T # sort rows according to y?
){
  # select and rename relevant columns for the plot   
  plot_df <- input_df %>% select(x_col=all_of(x_col),
                                        y_col=all_of(y_col),
                                        size_col=any_of(size_col),
                                        shape_col=any_of(shape_col),
                                        color_col=any_of(color_col),
                                        fill_col=any_of(fill_col),
                                        text_col=any_of(text_col),
                                        facet_col=any_of(facet_col))
  
  # creation of labels, abbreviating text 
  plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
  transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
  plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
  rm(transient)
  
  # order rows by x, remove eventual duplicates
  if(sort_y){
    plot_df<- plot_df %>% arrange(-x_col)
  }
  
  plot_df<-plot_df[!duplicated(plot_df %>% select(any_of(c("y_col","facet_col")))),]
  plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)
  
  lp <- ggplot(plot_df, aes(x_col, y_col)) +
    
    geom_segment(aes(x = 0, y = y_col, xend = x_col, yend = y_col))+
    theme_tufte(base_size = 16)+
    # scale_color_manual(values = "grey25", guide = "none")+
    theme(axis.title.y = element_blank(),axis.ticks.y= element_blank())+
    xlab(x_col)
  
  if("size_col" %in% colnames(plot_df)){
    lp <- lp + geom_point(aes(size=size_col)) + scale_size(name=size_col,range = size_vec)
  } else {
    lp <- lp + geom_point(size=2)      
  }
  
  if("shape_col" %in% colnames(plot_df)){
    lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
  }
  
  if("color_col" %in% colnames(plot_df)){
    lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) # drop = FALSE
  }
  
  if("fill_col" %in% colnames(plot_df)){
    lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
  }
  
  if("text_col" %in% colnames(plot_df)){
    lp <- lp + geom_text(aes(label=text_col, hjust=(-0.5)),fontface="italic", size=size_text, show.legend = F) +
      scale_x_continuous(expand = expansion(mult = c(.01, .1)))
  }
  
  if("facet_col" %in% colnames(plot_df)){
    lp <- lp + facet_col(vars(facet_col), scales = "free_y", space = "free") + 
      theme(strip.text = element_text(margin = margin(1,1,1,1), face="bold"))
  }
  
  return(lp)
  
}


### The function enrichment_b2b_lollipop creates a back to back enrichment plot (lollipop) starting from an enrichment dataframe with 2 analyses ----

enrichment_b2b_lollipop <- function(input_df, # input dataframe
                                    
                                    bb_col="input_name", # column used for the b2b (2 levels)
                                    x_col="combined_score", # column for x axis, necessary
                                    y_col="anno_name", # column for y axis, necessary
                                    
                                    size_col="overlap_size", # column for size, optional
                                    size_vec=c(1,5), # size vector c(min,max), optional
                                    
                                    shape_col="P<0.05", # column for shape, optional
                                    shape_vec=c(16,21), # shape vector, optional
                                    
                                    color_col="input_name", # column for color, optional
                                    color_vec=c("#16697A","#DB6400"), # color vector, optional
                                    
                                    fill_col="P<0.05", # column for fill color, optional
                                    fill_vec=c("grey25","white"), # fill vector, optional
                                    
                                    text_col="overlap_size", # column for fill color, optional
                                    
                                    facet_col="anno_class", # column for facet, optional
                                    
                                    char_max = 50, # maximum number of chars in text
                                    sort_y = T, # sort rows according to y?
                                    shift_mult = 1, # space left for middle labels (by default, if left to 1, equal to the space dedicated to the lollipop on each side)
                                    break_vec=seq(0,200,50) # specification of x breaks
){
  
  plot_df <- input_df %>% select(bb_col=all_of(bb_col),
                                        x_col=all_of(x_col),
                                        y_col=all_of(y_col),
                                        size_col=any_of(size_col),
                                        shape_col=any_of(shape_col),
                                        color_col=any_of(color_col),
                                        fill_col=any_of(fill_col),
                                        text_col=any_of(text_col),
                                        facet_col=any_of(facet_col))
  
  # creation of labels, abbreviating text 
  plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
  transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
  plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
  rm(transient)
  
  # creation of labels, abbreviating text 
  plot_df$x_col<- plot_df %>% pull(x_col) %>% str_trim() # Clean text
  transient <- plot_df$x_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
  plot_df$x_col<- ifelse(str_length(plot_df$x_col)>char_max,transient,plot_df$x_col) # replacement
  rm(transient)
  
  # factorize bb_col
  if(!is.factor(plot_df$bb_col)){
    plot_df$bb_col<-factor(plot_df$bb_col)
  }
  
  # sort y_col
  if(sort_y){
    plot_df<- plot_df %>% arrange(-x_col)
  }
  
  # remove eventual duplicates and factorize y_col 
  plot_df<-plot_df[!duplicated(plot_df %>% select(any_of(c("bb_col","y_col","facet_col")))),]
  plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)
  
  # define shift (if 1, roughly 1/3 of the plot)
  shift<-(max(plot_df$x_col)/2)*shift_mult
  
  f_left<- levels(plot_df$bb_col)[1]
  f_right<- levels(plot_df$bb_col)[2]
  
  lp<- ggplot(plot_df, aes(y = y_col))+
    geom_linerange(data = plot_df[plot_df$bb_col==f_left,],
                   aes(xmin = -shift, xmax = -shift-x_col))+
    geom_linerange(data = plot_df[plot_df$bb_col==f_right,],
                   aes(xmin = shift, xmax = shift+x_col))+
    geom_text(data = plot_df[plot_df$bb_col==f_left,],
              aes(y = y_col, x = 0, label = y_col),
              inherit.aes = F,size=3)+
    scale_x_continuous(limits =c((-shift-max(plot_df$x_col)),(shift+max(plot_df$x_col))),
                       breaks = c(rev(-break_vec)-shift, break_vec+shift),
                       labels = c(as.character(rev(break_vec)),as.character(break_vec)))+
    theme_tufte(base_size = 16)+
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top")+
    labs(x = x_col)
  
  
  if("size_col" %in% colnames(plot_df)){
    lp <- lp + 
      geom_point(data = plot_df[plot_df$bb_col==f_left,],
                 aes(x = -shift-x_col,size=size_col))+
      geom_point(data = plot_df[plot_df$bb_col==f_right,],
                 aes(x = shift+x_col,size=size_col))+
      scale_size(name=size_col,range = size_vec)
  } else {
    lp <- lp + 
      geom_point(data = plot_df[plot_df$bb_col==f_left,],
                 aes(x = -shift-x_col),size=2)+
      geom_point(data = plot_df[plot_df$bb_col==f_right,],
                 aes(x = shift+x_col),size=2)
  }
  
  if("text_col" %in% colnames(plot_df)){
    lp <- lp + geom_text(data = plot_df[plot_df$bb_col==f_left,],
                         aes(x = -shift-x_col, label=text_col, hjust=(1.1)),
                         fontface="italic", size=3, show_guide = F) +
      geom_text(data = plot_df[plot_df$bb_col==f_right,],
                aes(x = shift+x_col, label=text_col, hjust=(-0.1)),
                fontface="italic", size=3, show_guide = F) 
    scale_x_continuous(limits =c((-shift-max(plot_df$x_col)),(shift+max(plot_df$x_col))),
                       breaks = c(rev(-break_vec)-shift, break_vec+shift),
                       labels = c(as.character(rev(break_vec)),as.character(break_vec)),
                       expand = expansion(mult = c(.3, .3)))
  }
  
  if("shape_col" %in% colnames(plot_df)){
    lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
  }
  
  if("color_col" %in% colnames(plot_df)){
    lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) # drop = FALSE
  }
  
  if("fill_col" %in% colnames(plot_df)){
    lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
  }
  
  if("facet_col" %in% colnames(plot_df)){
    lp <- lp + facet_col(vars(facet_col), scales = "free_y", space = "free") + 
      theme(strip.text = element_text(margin = margin(1,1,1,1), face="bold"))
  }
  
  return(lp)
  
}



### The function deps_b2b_lollipop creates a back to back deps plot (lollipop) starting from an deps dataframe with 2 analyses ----

deps_b2b_lollipop <- function(input_df, # input dataframe
                              
                              bb_col="class", # column used for the b2b (2 levels)
                              x_col="N", # column for x axis, necessary
                              y_col="comp", # column for y axis, necessary
                              
                              size_col="N", # column for size, optional
                              size_vec=c(1,5), # size vector c(min,max), optional
                              
                              shape_col="class", # column for shape, optional
                              shape_vec=c(25,24), # shape vector, optional
                              
                              color_col="id", # column for color, optional
                              color_vec=c("grey60","grey40"), # color vector, optional
                              
                              fill_col="id", # column for fill color, optional
                              fill_vec=c("grey60","grey40"), # fill vector, optional
                              
                              text_col="N", # column for fill color, optional
                              
                              facet_col=c("Down-regulated", "Up-regulated"), # column for facet, optional
                              size_text=3,
                              char_max = 50, # maximum number of chars in text
                              sort_y = T, # sort rows according to y?
                              shift_mult = 1, # space left for middle labels (by default, if left to 1, equal to the space dedicated to the lollipop on each side)
                              break_vec=seq(0,1500,250), # specification of x breaks
                              position_dodge = 0.06
){
  bf<-"Arial"
  plot_df <- input_df %>% select(bb_col=all_of(bb_col),
                                 x_col=all_of(x_col),
                                 y_col=all_of(y_col),
                                 size_col=any_of(size_col),
                                 shape_col=any_of(shape_col),
                                 color_col=any_of(color_col),
                                 fill_col=any_of(fill_col),
                                 text_col=any_of(text_col),
                                 facet_col=any_of(facet_col))
  
  # creation of labels, abbreviating text 
  plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
  transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
  plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
  rm(transient)
  
  # factorize bb_col
  if(!is.factor(plot_df$bb_col)){
    plot_df$bb_col<-factor(plot_df$bb_col)
  }
  
  # sort y_col
  if(sort_y){
    plot_df<- plot_df %>% arrange(-x_col)
  }
  
  # remove eventual duplicates and factorize y_col 
  plot_df<-plot_df[!duplicated(plot_df %>% select(any_of(c("bb_col","y_col","facet_col")))),]
  plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)
  
  # define shift (if 1, roughly 1/3 of the plot)
  shift<-(max(plot_df$x_col)/2)*shift_mult
  
  if(uniqueN(plot_df$bb_col) != 1){
    f_left<- levels(plot_df$bb_col)[1]
    f_right<- levels(plot_df$bb_col)[2]
    
    lp<- ggplot(plot_df, aes(y = y_col))+
      geom_linerange(data = plot_df[plot_df$bb_col==f_left,],
                     aes(xmin = -shift, xmax = -shift-x_col, size=0.75))+
      geom_linerange(data = plot_df[plot_df$bb_col==f_right,],
                     aes(xmin = shift, xmax = shift+x_col, size=0.75))+
      geom_text(data = plot_df[plot_df$bb_col==f_left,],
                aes(y = y_col, x = 0, label = y_col),
                inherit.aes = F,size=size_text)+
      geom_text(data = plot_df[plot_df$bb_col==f_right,],
                aes(y = y_col, x = 0, label = y_col),
                inherit.aes = F,size=size_text)+
      scale_x_continuous(limits =c((-shift-max(break_vec)),(shift+max(break_vec))),
                         breaks = c(rev(-break_vec)-shift, break_vec+shift),
                         labels = c(as.character(rev(break_vec)),as.character(break_vec)))+
      theme_tufte(base_size = 16)+
      theme(axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none")+
      labs(x = "N DEPs")
    
    
    if("size_col" %in% colnames(plot_df)){
      lp <- lp + 
        geom_point(data = plot_df[plot_df$bb_col==f_left,],
                   aes(x = -shift-x_col, size = 5), position = position_nudge(y = -position_dodge))+
        geom_point(data = plot_df[plot_df$bb_col==f_right,],
                   aes(x = shift+x_col, size = 5), position = position_nudge(y = position_dodge))+
        scale_size(name=size_col,range = size_vec)
    } else {
      lp <- lp + 
        geom_point(data = plot_df[plot_df$bb_col==f_left,],
                   aes(x = -shift-x_col),size=2)+
        geom_point(data = plot_df[plot_df$bb_col==f_right,],
                   aes(x = shift+x_col),size=2)
    }
    
    
    if("text_col" %in% colnames(plot_df)){
      lp <- lp + geom_text(data = plot_df[plot_df$bb_col==f_left,],
                           aes(x = -(shift*13/10)-x_col, label=text_col, vjust=0),
                           fontface="italic", size=size_text, show.legend = F) +
        geom_text(data = plot_df[plot_df$bb_col==f_right,],
                  aes(x = (shift*13/10)+x_col, label=text_col, vjust=0),
                  fontface="italic", size=size_text, show.legend = F)
    }
    
    if("shape_col" %in% colnames(plot_df)){
      lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
    }
    
    if("color_col" %in% colnames(plot_df)){
      lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col,values= color_vec) # drop = FALSE
    }
    
    if("fill_col" %in% colnames(plot_df)){
      lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
    }
    
    if(uniqueN(plot_df$y_col) == length(unique(plot_df$y_col))){
      lp <- lp + geom_text(aes(x=shift, label=c("Up-regulated", rep("",nrow(plot_df)-1)), hjust=(-0.05), vjust=(-2), fontface=2, color="black"),size=size_text)
      lp <- lp + geom_text(aes(x=-shift, label=c("Down-regulated", rep("",nrow(plot_df)-1)), hjust=(1.1), vjust=(-2), fontface=2, color="black"),size=size_text)
    } else{
      lp <- lp + geom_text(aes(x=shift, label=c("Up-regulated", rep("",nrow(plot_df)-1)), hjust=(-0.05), vjust=(-2), fontface=2, color="black"),size=size_text, nudge_y = 2.5)
      lp <- lp + geom_text(aes(x=-shift, label=c("Down-regulated", rep("",nrow(plot_df)-1)), hjust=(1.1), vjust=(-2), fontface=2, color="black"),size=size_text, nudge_y = 2.5)
    }
  } else{ # Plot for interactomics
    f_right<- levels(plot_df$bb_col)[1]
    
    lp<- ggplot(plot_df, aes(y = y_col))+
      geom_linerange(data = plot_df[plot_df$bb_col==f_right,],
                     aes(xmin = shift, xmax = shift+x_col, size=0.75))+
      geom_text(data = plot_df[plot_df$bb_col==f_right,],
                aes(y = y_col, x = 0, label = y_col),
                inherit.aes = F,size=size_text)+
      scale_x_continuous(limits =c(-300,(shift+max(break_vec))),
                         breaks = c(break_vec+shift),
                         labels = c(as.character(break_vec)))+
      theme_tufte(base_size = 16)+
      theme(axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none")+
      labs(x = "N DEPs")
    
    
    if("size_col" %in% colnames(plot_df)){
      lp <- lp + 
        geom_point(data = plot_df[plot_df$bb_col==f_right,],
                   aes(x = shift+x_col, size = 5), position = position_nudge(y = position_dodge))+
        scale_size(name=size_col,range = size_vec)
    } else {
      lp <- lp + 
        geom_point(data = plot_df[plot_df$bb_col==f_right,],
                   aes(x = shift+x_col),size=2)
    }
    
    
    if("text_col" %in% colnames(plot_df)){
      lp <- lp + 
        geom_text(data = plot_df[plot_df$bb_col==f_right,],
                  aes(x = (shift*13/10)+x_col, label=text_col, vjust=0),
                  fontface="italic", size=size_text, show.legend = F)
    }
    
    if("shape_col" %in% colnames(plot_df)){
      lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec) # drop = FALSE
    }
    
    if("color_col" %in% colnames(plot_df)){
      lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col,values= color_vec) # drop = FALSE
    }
    
    if("fill_col" %in% colnames(plot_df)){
      lp <- lp + aes(fill=fill_col) + scale_fill_manual(name=fill_col, values= fill_vec) # drop = FALSE
    }
    
    if(nrow(plot_df) == length(unique(plot_df$y_col))){
      lp <- lp + geom_text(aes(x=shift, label=c("Up-regulated", rep("",nrow(plot_df)-1)), hjust=(-0.05), vjust=(-2), fontface=2, color="black"),size=size_text)
    } else{
      lp <- lp + geom_text(aes(x=shift, label=c("Up-regulated", rep("",nrow(plot_df)-1)), hjust=(-0.05), vjust=(-2), fontface=2, color="black"),size=size_text, nudge_y = 2.5)
    }
  }
  
  return(lp)
  
}




### The function enrichment_dotmatrix creates am enrichment matrix dotplot, starting from an enrichment dataframe with multiple analyses ----

enrichment_dotmatrix <- function(input_df, # input dataframe
                                 x_col="input_name", # column for x axis, necessary
                                 y_col="anno_name", # column for y axis, necessary
                                 
                                 size_col="odds_ratio", # column for size, optional
                                 size_vec=c(1,5), # size vector c(min,max), optional
                                 
                                 shape_col="P<0.05", # column for shape, optional
                                 shape_vec=c(16,21), # shape vector, optional
                                 
                                 color_col="P<0.05", # column for color, optional
                                 color_vec=c("grey25","grey50"), # color vector, optional
                                 
                                 fill_col="P<0.05", # column for fill color, optional
                                 fill_vec=c("grey25","white"), # fill vector, optional
                                 
                                 facet_col="anno_class", # column for facet, optional
                                 
                                 char_max = 50, # maximum number of chars in text (y_axis)
                                 sort_y = F, # sort x rows? (in development)
                                 sort_x = F # sort y rows? according to size_col (in development)
){
  input_df$color<-ifelse(input_df$Significant==TRUE, color_vec[input_df$input_name], "white")
  
  # select and rename relevant columns for the plot   
  plot_df <- input_df %>% select(x_col=all_of(x_col),
                                        y_col=all_of(y_col),
                                        size_col=any_of(size_col),
                                        shape_col=any_of(shape_col),
                                        color_col=any_of(color_col),
                                        fill_col=any_of(fill_col),
                                        facet_col=any_of(facet_col),
                                        color=any_of("color"))
  
  # creation of labels, abbreviating text 
  plot_df$y_col<- plot_df %>% pull(y_col) %>% str_trim() # Clean text
  transient <- plot_df$y_col %>% str_sub(1,(char_max-2)) %>% str_c("..") # abbreviated option
  plot_df$y_col<- ifelse(str_length(plot_df$y_col)>char_max,transient,plot_df$y_col) # replacement
  rm(transient)
  
  # order rows by x, remove eventual duplicates
  if(sort_y | sort_x){
    plot_df<- plot_df %>% arrange(-size_col)
  }
  
  plot_df<-plot_df[!duplicated(plot_df %>% select(any_of(c("x_col","y_col","facet_col")))),]
  
  if(!(is.factor(plot_df$x_col))){
    plot_df$x_col<-factor(plot_df$x_col,levels=(unique(plot_df$x_col)))# factorization (x_col)
  }
  if(!(is.factor(plot_df$y_col))){
    plot_df$y_col<-factor(plot_df$y_col,levels=rev(unique(plot_df$y_col))) # factorization (y_col)
  }
  
  
  lp <- ggplot(plot_df, aes(x_col, y_col)) +
    theme_tufte(base_size = 16)+
    theme(axis.title = element_blank(),axis.text.x = element_text(angle = 30, hjust = 1))  
  
  if("size_col" %in% colnames(plot_df)){
    lp <- lp + geom_point(aes(size=size_col)) + scale_size(name=size_col,range = size_vec)
  } else {
    lp <- lp + geom_point(size=2)      
  }
  
  if("shape_col" %in% colnames(plot_df)){
    lp <- lp + aes(shape=shape_col) + scale_shape_manual(name=shape_col, values= shape_vec, labels=if("UP" %in% plot_df$shape_col){c("TRUE","FALSE","")}else{unique(plot_df$shape_col)}) # drop = FALSE
  }
  
  if("color_col" %in% colnames(plot_df)){
    lp <- lp + aes(colour=color_col) + scale_color_manual(name=color_col, values= color_vec) + guides(color="none") # drop = FALSE
  }
  
  if("fill_col" %in% colnames(plot_df)){
    cc<-plot_df$color
    names(cc)<-interaction(plot_df$x_col,plot_df$y_col,sep="-",lex.order=TRUE)
    lp <- lp + aes(fill=interaction(x_col,y_col,sep="-",lex.order=TRUE)) + scale_fill_manual(name="Significant", values= cc) + guides(fill="none") # drop = FALSE
  }
  
  if("facet_col" %in% colnames(plot_df)){
    lp <- lp + facet_col(vars(facet_col), scales = "free_y",  space = "free") +
      theme(strip.text = element_text( margin = margin(1,1,1,1)))
  }
  
  return(lp)
  
}



### The function enrichment_geneterm_network creates a gene-term network starting from an enrichment dataframe ----

enrichment_geneterm_network <- function(input_df){
  
  gene_col <- "overlap_ids"
  term_col <- "anno_name"
  
  col_vec<- c("#B2473E","black")
  
  plot_df <- input_df %>% select(term_col=all_of(term_col),
                                        gene_col=all_of(gene_col))
  
  edges_df <- plot_df  %>% separate_rows(gene_col,sep=";")
  edges_df$gene_col <- str_to_title(edges_df$gene_col) %>% str_trim()
  
  plot_df$term_col <- str_split(plot_df$term_col,"\\(",simplify=T)[,1] %>% str_trim()
  
  
  terms_df<- tibble(name=unique(edges_df$term_col),
                    class="term")
  genes_df<- tibble(name=unique(edges_df$gene_col),
                    class="gene")
  nodes_df<-bind_rows(terms_df,genes_df)
  
  net_net <- fortify(as.edgedf(as.data.frame(edges_df)), nodes_df)
  net_net <- subset(net_net, from_id!=to_id)
  g <- graph_from_data_frame(net_net, directed=FALSE, vertices = nodes_df)
  l<-ggraph(g,layout="kk") +
    geom_edge_link2(aes(edge_colour = "grey50",alpha=0.2)) +
    ggtitle("Network",subtitle = paste0("Layout: kk"))+
    scale_edge_width(range = c(0.1,1)) +
    scale_edge_colour_manual(values = "grey50") +
    geom_node_point(aes(fill=class), shape = 21,size=2.5,pch='.') +
    geom_node_label(aes(filter = class=="term", label = name, size=2), repel=TRUE) +
    geom_node_text(aes(filter = class=="gene", label = name, size=2), repel=TRUE,check_overlap = TRUE) +
    scale_fill_manual(values = col_vec) +
    scale_color_brewer(palette ="Dark2") +
    theme_graph() +
    theme(legend.position = "bottom") +
    guides(text=F)
  l
}

### Function to extract coverage data from bedtools coverage with option d (coverage for each nucleotide) ----


resume_bedcoverage_d<-function(dir,
                               tailcut=4) {
  
  list_input<-list.files(path=dir, pattern="bed", full.names=F)
  
  if(exists("o_tab")){
    rm(o_tab)
  }
  
  for(i in (list_input)){
    i_name<-substr(i,1,(nchar(i)-tailcut))
    i_tab<-read.table(paste0(dir,i),header=FALSE,sep="\t")
    if(!exists("o_tab")){
      o_tab<-data.frame(name=i_tab[,8])
      rownames(o_tab)<-paste(i_tab[,4],i_tab[,7],sep="_") # correct end (-1)
      colnames(o_tab)<-c(i_name)
    } else {
      o_tab[,i_name]<-i_tab[,8]
    }
  }
  return(o_tab)
}


### Function to use edgeR function glmQLF generalized linear models (quasi likelihood function) ----
edgeRglmQLF<-function(mat=edge_f, # object of class DGEGLM
                      contro, # comparison, create with makeContrasts
                      cpm_mat=edge_n, #used to calculate average signal
                      label="", # include in label column
                      sig_thr=0.5, # signal threshold
                      sig_col="CPM", # signal column
                      fc_thr=0.5, # fold change threshold
                      pval_thr=0.05, # p-value threshold
                      pval_col="p_val", # p-value column
                      names=FALSE)
{ # include label in column names
  degs<-glmQLFTest(edge_f,contrast=contro)$table[,-3]
  colnames(degs)<-c("log2_FC","log2_CPM","p_val")
  a_levels<-rownames(contro)[which(contro!=0)]
  a_samples<-which(cpm_mat$samples$group%in%a_levels)
  cpm_sele<-cpm(cpm_mat,log=T)[,a_samples]
  degs$log2_CPM<-apply(cpm_sele,1,function(x) mean(x))
  degs$p_adj<-p.adjust(degs$p_val, method ="BH")
  degs$CPM<-round(2^degs$log2_CPM,2)
  degs$class<-"="
  degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC>=fc_thr & degs[,pval_col]<=pval_thr),"class"]<-"+"
  degs[which(degs[,sig_col]>=sig_thr & degs$log2_FC<=(-fc_thr) & degs[,pval_col]<=pval_thr),"class"]<-"-"
  degs$class<-as.factor(degs$class)
  degs$comp<-label
  degs$id<-rownames(degs)
  degs<-degs[,c("id","comp","log2_FC","CPM","p_val","p_adj","class")]
  if(names=="TRUE"){
    newnames<-paste(label,colnames(degs),sep="_")
    colnames(degs)<-newnames
  }
  return(degs)
}

### Function to extract and save only the different part of two string ----
mf <- function(x,sep){
  xsplit = strsplit(x,split = sep)
  xdfm <- as.data.frame(do.call(rbind,xsplit))
  res <- list()
  for (i in 1:ncol(xdfm)){
    if (!all(xdfm[,i] == xdfm[1,i])){
      res[[length(res)+1]] <- as.character(xdfm[,i])
    }
  }
  res <- as.data.frame(do.call(rbind,res))
  res <- apply(res,2,function(x) paste(x,collapse="_"))
  return(res)
}

### Function for the enrichment of the protein in the dataset. It use the function enrichment_enrichr ----
enrichRfnc<-function(in_df, pval_fdr_enrich, pval_enrich_thr, overlap_size_enrich_thr, dbs=NULL, with_background=FALSE){
  setEnrichrSite("Enrichr") 
  DEGs_lists<-NULL
  background_lists<-NULL
  
  for(comp_c in unique(in_df$comp)){
    
    DEGs_lists[[paste0(comp_c,"_all")]]<-in_df %>% filter(comp==comp_c,class!="=") %>% pull(id) %>% unique() %>% sort()
    DEGs_lists[[paste0(comp_c,"_up")]]<-in_df %>% filter(comp==comp_c,class=="+") %>% pull(id) %>% unique() %>% sort()
    DEGs_lists[[paste0(comp_c,"_down")]]<-in_df %>% filter(comp==comp_c,class=="-") %>% pull(id) %>% unique() %>% sort()
    
    if(with_background){
      background_lists[[paste0(comp_c,"_all")]]<-in_df %>% filter(comp==comp_c) %>% pull(id) %>% unique() %>% sort()
      background_lists[[paste0(comp_c,"_up")]]<-in_df %>% filter(comp==comp_c) %>% pull(id) %>% unique() %>% sort()
      background_lists[[paste0(comp_c,"_down")]]<-in_df %>% filter(comp==comp_c) %>% pull(id) %>% unique() %>% sort()
    }
  }
  ncores <- min(detectCores()-2, length(unique(in_df$comp)))
  enr_df<-NULL
  back_frg<-NULL
  
  if(.Platform$OS.type == "unix") {
    enrfcn <- function(a) {
      # source("R/functions.R")
      frg<-DEGs_lists[[a]]
      if(with_background){back_frg <- background_lists[[a]]}
      if(length(frg)>0){
        enrichment_enrichr(frg, input_name=a, dbs_vec = dbs, back_frg = back_frg)
      }
    }
    enr_df<-do.call(rbind, mclapply(names(DEGs_lists), enrfcn, mc.cores = ncores))
  } else {
    cluster_ext <- makeCluster(ncores, type = "SOCK")
    registerDoParallel(cl = cluster_ext)
    
    enr_df<-foreach(a = names(DEGs_lists),.combine='rbind',.packages = c("dplyr","enrichR","tidyr","stringr","proTN")) %dopar% {
      # source("R/functions.R")
      frg<-DEGs_lists[[a]]
      if(with_background){back_frg <- background_lists[[a]]}
      if(length(frg)>0){
        enrichment_enrichr(frg, input_name=a, dbs_vec = dbs, back_frg = back_frg)
      }
    }
    stopCluster(cluster_ext)
  }
  
  if(any(is.na(enr_df))){
    enr_df<-NULL
    for(a in names(DEGs_lists)){
      frg<-DEGs_lists[[a]]
      if(with_background){back_frg <- background_lists[[a]]}
      if(length(frg)>0){
        enr_df<-rbind(enr_df,enrichment_enrichr(input_vec = frg, input_name=a, dbs_vec = dbs, back_frg = back_frg))
      }
    }
  }
  enr_df<-na.omit(enr_df)
  enr_df<-
    transform(enr_df, 
              overlap_size = as.numeric(overlap_size), 
              p_value = as.numeric(p_value), 
              fdr = as.numeric(fdr), 
              odds_ratio = as.numeric(odds_ratio), 
              combined_score = as.numeric(combined_score), 
              input_size = as.numeric(input_size), 
              anno_size = as.numeric(anno_size), 
              overlap_input_ratio = as.numeric(overlap_input_ratio),
              overlap_anno_ratio = as.numeric(overlap_anno_ratio))
  enr_df$"Significant"<-ifelse(if(pval_fdr_enrich=="p_adj"){
    (enr_df$fdr<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)
    }else if(pval_fdr_enrich=="p_val"){
      (enr_df$p_value<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)
      },"TRUE","FALSE") %>% factor(levels=c("TRUE","FALSE"))
  enr_df$log2_OR<-log2(enr_df$odds_ratio)
  return(enr_df)
}

### Function for the enrichment of the protein of an universe. It use the function enrichment_enrichr ----
enrichRfnc_universe<-function(in_df, pval_fdr_enrich, pval_enrich_thr, overlap_size_enrich_thr, dbs=NULL){
  enr_df <- tryCatch({
    DEGs_lists<-NULL
    DEGs_lists[["Universe_all"]]<-in_df
    
    enr_df<-enrichment_enrichr(DEGs_lists[["Universe_all"]], input_name="Universe_all", dbs_vec = dbs)
    
    enr_df$"Significant"<-ifelse(if(pval_fdr_enrich){(enr_df$fdr<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)}else{(enr_df$p_value<pval_enrich_thr) & (enr_df$overlap_size>=overlap_size_enrich_thr)},"TRUE","FALSE") %>% factor(levels=c("TRUE","FALSE"))
    enr_df$log2_OR<-log2(enr_df$odds_ratio)
    enr_df <- enr_df[-which(enr_df$Significant == "FALSE"),]
    enr_df
  },
  error=function(cond){
    print("\n ERROR: An error occur when connect to EnrichR. \n ")
    NULL
  })
  return(enr_df)
}


### Update limma function with data.table ----
limmafnc_dt <- function(type = "PROT", c_anno, dat_gene, psm_count_table, formule_contrast, expr_avgse_df, signal_thr, fc_thr, pval_thr, pval_fdr, interactomics = FALSE) {
  # make design table
  design <- model.matrix(~0 + c_anno$condition)
  colnames(design) <- levels(as.factor(c_anno$condition))
  rownames(design) <- c_anno$sample
  
  gene_matrix <- as.data.frame(dat_gene[, rownames(design), with = FALSE], row.names = as.character(unlist(dat_gene[,1])))
  gene_matrix <- as.matrix(gene_matrix)
  
  # DEqMS part analysis
  setnames(psm_count_table, c("gene", "count"))
  setkey(psm_count_table, gene)
  
  filt_contro_list <- Filter(function(i) {
    all(str_extract_all(formule_contrast[i], "\\w+")[[1]] %in% c(colnames(design), as.character(0:10)))
  }, seq_along(formule_contrast))
  
  formule_contrast <- formule_contrast[unlist(filt_contro_list)]
  if (length(formule_contrast) == 0) stop("Error: No valid contrast design given.")
  
  contrast <- makeContrasts(contrasts = formule_contrast, levels = design)
  colnames(contrast) <- names(formule_contrast)
  
  # limma part analysis
  fit1 <- lmFit(gene_matrix, design)
  fit2 <- eBayes(contrasts.fit(fit1, contrasts = contrast))
  
  if (type == "PROT") {
    fit2$count <- psm_count_table[match(gene, rownames(fit2$coefficients)), count]
    fit3 <- spectraCounteBayes(fit2)
  } else {
    fit3 <- fit2
  }
  
  # extract results
  degs_w_df <- data.table()
  degs_l_df <- data.table()
  
  signal_col <- "log2_expr"
  fc_col <- "log2_FC"
  pval_col <- pval_fdr
  
  for (comp in (colnames(fit3$coefficients))) {
    if(type == "PROT"){
      degs_u<-tryCatch({
        outputResult(fit3,coef_col = comp) %>% 
          rename("log2_FC"="logFC",
                        "log2_expr"="AveExpr",
                        "p_val"="sca.P.Value",
                        "p_adj"="sca.adj.pval",
                        "id"="gene") %>% arrange(id)
      }, 
      error=function(cond){
        colnames(fit3$sca.t) <- colnames(fit3$coefficients)
        colnames(fit3$sca.p) <- colnames(fit3$coefficients)
        outputResult(fit3,coef_col = comp) %>% 
          rename("log2_FC"="logFC",
                        "log2_expr"="AveExpr",
                        "p_val"="sca.P.Value",
                        "p_adj"="sca.adj.pval",
                        "id"="gene") %>% arrange(id)
      })
    }else if(type == "PEP"){
      degs_u<-topTable(fit3, coef = comp, number = Inf) %>% 
        rename("log2_FC"="logFC",
                      "log2_expr"="AveExpr",
                      "p_val"="P.Value",
                      "p_adj"="adj.P.Val") %>% mutate("id"=rownames(topTable(fit3, coef = comp, number = Inf))) %>% select(-c(t,B))
    }
    
    degs_u <- as.data.table(degs_u)
    if(!interactomics){
      degs_u[, `:=`(log2_expr = 2^log2_FC, FC = 2^log2_FC, comp = comp, class = "=")]
      degs_u[log2_FC >= fc_thr & get(pval_col) <= pval_thr, class := "+"]
      degs_u[log2_FC <= -fc_thr & get(pval_col) <= pval_thr, class := "-"]
    } else{
      degs_u[, `:=`(log2_expr = 2^log2_FC, FC = 2^log2_FC, comp = comp, class = "=")]
      degs_u[log2_FC >= fc_thr & get(pval_col) <= pval_thr, class := "+"]
    }
    
    degs_l_df <- rbindlist(list(degs_l_df, degs_u[, .(id, comp, class, log2_FC, FC, p_val, p_adj)]))
    degs_add <- degs_u[, .(id, class, log2_FC, FC, p_val, p_adj)]
    setnames(degs_add, names(degs_add)[-1], paste(comp, names(degs_add)[-1], sep = "_"))
    if(nrow(degs_w_df)==0){
      degs_w_df <- degs_add
    } else{
      degs_w_df <- merge(degs_w_df, degs_add, by = "id", all = TRUE)
    }
  }
  
  return(list(degs_l_df = degs_l_df, degs_w_df = degs_w_df))
}
