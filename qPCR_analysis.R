library(tidyverse)
library(ggpubr)

setwd("/home/zaczou/projects/CRISPR/qPCR")

#-----------------read data------------------#

data_1 <- read_tsv("data_1.txt", skip = 1, col_names = T)
panel_1 <- read_tsv("panel_1.1.dict.txt", col_names = T)

data_2.1 <- read_tsv("data_2.1.txt", skip = 1, col_names = T)
panel_2.1 <- read_tsv("panel_2.1.dict.txt", col_names = T)

data_2.2 <- read_tsv("data_2.2.txt", skip = 1, col_names = T)
panel_2.2 <- read_tsv("panel_2.2.dict.txt", col_names = T)

#-----------------------preprocessing-----------------#

exp_1 <- panel_1 %>%
  dplyr::left_join(data_1[,c("Pos","Cp")], by = c("position"="Pos")) %>%
  dplyr::mutate(position=NULL, group="M0") %>%
  na.omit()

exp_2.1 <- panel_2.1 %>%
  dplyr::left_join(data_2.1[,c("Pos","Cp")], by = c("position"="Pos"))
exp_2.2 <- panel_2.2 %>%
  dplyr::left_join(data_2.2[,c("Pos","Cp")], by = c("position"="Pos"))
exp_2 <- bind_rows(exp_2.1, exp_2.2) %>%
  dplyr::mutate(position=NULL, group="TAM") %>%
  na.omit()
exp <- bind_rows(exp_1, exp_2) %>%
  unite(id, c(primer, sample, group), remove = F)
for(ii in 1:nrow(exp)){if(exp$sample[ii]=="s11"){exp$group[ii] = "THP-1"}}
for(ii in 1:nrow(exp)){if(exp$sample[ii]=="s7"|exp$sample[ii]=="s10"){exp$group[ii] = "Mono_TAM"}}
for(ii in 1:nrow(exp)){if(exp$group[ii]=="M0"&exp$sample[ii]=="s6"){exp$group[ii] = "Monocyte"}}
exp <- exp %>% na.omit()

saveRDS(exp, "exp.rds")
exp <- readRDS("exp.rds")

#---------------------data laundry-----------------#

exp.lim <- exp %>% 
  dplyr::group_by(id) %>%
  dplyr::summarise(lower = quantile(Cp,0.5)-IQR(Cp),
                   upper = quantile(Cp,0.5)+IQR(Cp),
                   id = id,
                   n = n()) %>% unique()

sexp <- exp
for(ii in 1:nrow(exp))
{
  tmp.id <- exp.lim$id==exp$id[ii]
  tmp.lo <- exp.lim$lower[tmp.id]
  tmp.hi <- exp.lim$upper[tmp.id]
  tmp.cp <- exp$Cp[ii]
  if(exp.lim$n[tmp.id]==3&(tmp.cp>tmp.hi|tmp.cp<tmp.lo))
  {sexp$Cp[ii] <- NA}
}
sexp <- na.omit(sexp)

#check
check <- function(tb){
  tmp <- tb %>% group_by(id) %>% tally()
  table(tmp$n)}

saveRDS(sexp, "sexp.rds")
sexp <- readRDS("sexp.rds")

#-----------------------delta----------------------#

sexp <- sexp %>% unite(id_sg, c(sample, group), remove = F)

count <- sexp %>%
  dplyr::filter(primer == "p21"|primer == "p22") %>%
  dplyr::group_by(id_sg) %>%
  dplyr::summarize(mCp = mean(Cp))
dexp <- sexp[sexp$primer!="p21"&sexp$primer!="p22",]
for(ii in 1:nrow(dexp))
{dexp[ii, "Cp"] = dexp[ii, "Cp"] - count[count$id_sg == unlist(dexp[ii, "id_sg"]), "mCp"]}

dexp <- dexp %>% dplyr::mutate(id_sg=NULL)

saveRDS(dexp, "dexp.rds")
dexp <- readRDS("dexp.rds")

#-------------------------delta delta----------------------#

ctrl.thp <- dexp %>%
  dplyr::filter(group == "THP-1") %>%
  dplyr::group_by(primer) %>%
  dplyr::summarize(mCp = mean(Cp))
ctrl.mnc <- dexp %>%
  dplyr::filter(group == "Monocyte") %>%
  dplyr::group_by(primer) %>%
  dplyr::summarize(mCp = mean(Cp))

ddexp <- dexp
for(ii in 1:nrow(ddexp))
{
  if(dexp[ii, "group"]=="THP-1"|dexp[ii, "group"]=="Monocyte"){ddexp[ii, "Cp"] = 0}
  else if(dexp[ii, "group"]=="Mono_TAM")
  {
    if(length(table(ctrl.mnc$primer == unlist(dexp[ii, "primer"])))==1){ddexp[ii, "Cp"] = NA}
    else
    {ddexp[ii, "Cp"] = dexp[ii, "Cp"] - ctrl.mnc[ctrl.mnc$primer == unlist(dexp[ii, "primer"]), "mCp"]}
  }
  else
  {ddexp[ii, "Cp"] = dexp[ii, "Cp"] - ctrl.thp[ctrl.thp$primer == unlist(dexp[ii, "primer"]), "mCp"]}
}
ddexp <- na.omit(ddexp)

saveRDS(ddexp, "ddexp.rds")
ddexp <- readRDS("ddexp.rds")

#---------------------fold change----------------------#

fexp <- ddexp %>%
  dplyr::mutate(FC = 2^(-Cp))

saveRDS(fexp, "fexp.rds")
fexp <- readRDS("fexp.rds")

#---------------------annotation-----------------------------#

res <- fexp %>%
  dplyr::left_join(marker, by = "primer") %>% 
  dplyr::left_join(cond, by = c("sample","group"))

marker <- read_tsv("marker.txt", col_names = T)

cond1 <- read_tsv("cond_1.txt", col_names = T) %>% dplyr::mutate(group="M0")
cond2 <- read_tsv("cond_2.txt", col_names = T) %>% dplyr::mutate(group="TAM")
cond3 <- tibble(sample = c("s7","s10","s6"),
                condition = c("Treated","Treated","Control"),
                group = c("Mono_TAM","Mono_TAM","Monocyte"))
cond4 <- tibble(sample = "s11", condition = "Control", group = "THP-1")
cond <- bind_rows(cond1, cond2, cond3, cond4)

mseq <- read_tsv("marker-seq.txt", col_names = T)
cseq <- read_tsv("cond-seq.txt", col_names = T)

res$condition <- factor(res$condition, levels = c("Control","Treated",cseq$condition))
res$group <- factor(res$group, levels = c("THP-1","M0","TAM","Monocyte","Mono_TAM"))
res$gene <- factor(res$gene, levels = mseq$gene)

saveRDS(res, "res.rds")
res <- readRDS("res.rds")

#--------------------------plot----------------------------#

# plot1 each gene

plot1sg <- function(cell, pr.gene)
{
  name <- paste(cell,"_",pr.gene,".png", sep = '')
  print(name)
  tmp <- res %>%
    dplyr::filter((
      (group==cell|group=="THP-1")&(primer==marker$primer[marker$gene==pr.gene])
    )) %>%
    na.omit() %>%
    group_by(condition) %>%
    summarize(ave_FC = mean(FC),
              se = sd(FC)/sqrt(n()),
              condition = condition,
              group = group) %>% unique()
  tmp %>% ggplot(aes(condition, ave_FC)) +
    geom_col(aes(fill = group), width = 0.8) +
    geom_errorbar(data = tmp[-1,], 
                  aes(x = condition, y = ave_FC, 
                      ymin=ave_FC-se, 
                      ymax=ave_FC+se),
                  width=.5,
                  position=position_dodge(.9), 
                  color="black") +
    scale_y_continuous(sec.axis = sec_axis(~., breaks = c(0,1.0))) +
    ylab("Relative RNA Expression") + xlab("Condition") +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    ggtitle(paste(cell, pr.gene, sep = ' ')) +
    theme(legend.title = element_blank()) +
    guides(x = guide_axis(angle = 30))
  ggsave(width = 5, height = 4, device = "png", name, dpi = 300, path = "./results/group1/single/")
}

plot1pl <- function(pr.gene)
{
  name <- paste("comp_",pr.gene,".png", sep = '')
  print(name)
  tmp <- res %>%
    dplyr::filter((
      (condition %in% cond1$condition|group=="THP-1")&(primer==marker$primer[marker$gene==pr.gene])
    )) %>%
    na.omit() %>%
    unite(tag, c(condition, group), remove = F) %>%
    group_by(tag) %>%
    summarize(ave_FC = mean(FC),
              se = sd(FC)/sqrt(n()),
              condition = condition,
              group = group) %>% unique()
  tmp %>% ggplot(aes(condition, ave_FC, fill = group)) +
    geom_bar(position = position_dodge2(width = 1.2), stat = "identity", width = 0.8) +
    geom_errorbar(data = filter(tmp,condition!="Control"), 
                  aes(x = condition, y = ave_FC, 
                      ymin=ave_FC-se, 
                      ymax=ave_FC+se),
                  width=.5,
                  position=position_dodge(0.8), 
                  color="black") +
    scale_y_continuous(sec.axis = sec_axis(~., breaks = c(0,1.0))) +
    ylab("Relative RNA Expression") + xlab("Condition") +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    ggtitle(pr.gene) +
    theme(legend.title = element_blank()) +
    guides(x = guide_axis(angle = 30))
  ggsave(width = 4, height = 4, device = "png", name, dpi = 300, path = "./results/group1/compare")
}

# plot2 each condition

clim <- function(lims)
{
  lim1 <- lims[1]
  lim2 <- lims[2]
  lim2 <- lim1+1.25*(lim2-lim1)
  return(c(lim1,lim2))
}

plot2 <- function(con)
{
  name <- paste("(",con,").png", sep = '')
  print(name)
   tmp <- res %>%
    dplyr::filter((
      (condition==con)|
      (group=="THP-1"))
      &(!primer %in% c("p20"))
      ) %>%
    na.omit()
  tmp.er <- tmp %>%
    unite(tag, c(group, gene), remove=FALSE) %>%
    group_by(tag) %>%
    summarize(ave_FC = mean(FC),
              se = sd(FC)/sqrt(n()),
              gene = gene,
              group = group) %>% unique()
  tmp %>% ggplot(aes(group, FC)) +
    geom_col(data=tmp.er, aes(x=group, y=ave_FC, fill=group), width = 0.8) +
    geom_errorbar(data=filter(tmp.er,group!="THP-1"), aes(
      x=group, y=ave_FC, 
      ymin=ave_FC-se, ymax=ave_FC+se
      ), width=.5, position=position_dodge(.9)) +
    stat_compare_means(label = "p.signif", method = "anova", size = 3, label.y.npc = rep(0.9,nrow(tmp.er)), label.x = 1.7) +
    scale_y_continuous(limits = clim, sec.axis = sec_axis(~., breaks = c(0,1.0))) +
    ylab("Relative RNA Expression") + xlab("Target Gene") +
    facet_wrap(vars(gene), nrow = 3, scales = "free_y") +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    ggtitle(con) +
    theme(axis.text.x=element_blank(), 
          axis.text=element_text(size=6), 
          strip.text=element_text(size=6),
          legend.title = element_blank()
          )
  ggsave(width = 8, height = 4, units = "in", device = "png", name, dpi = 300, path = "./results/group2/linear/")
}

# plot3 monocyte

plot3 <- function()
{
  name <- paste("Monocyte.png", sep = '')
  print(name)
  tmp <- res %>%
    dplyr::filter((
      (group=="Monocyte")|
      (group=="Mono_TAM"))
      &(!primer %in% c("p20"))
    ) %>%
    na.omit()
  tmp.er <- tmp %>%
    unite(tag, c(group, gene), remove=FALSE) %>%
    group_by(tag) %>%
    summarize(ave_FC = mean(FC),
              se = sd(FC)/sqrt(n()),
              gene = gene,
              group = group) %>% unique()
  tmp %>% ggplot(aes(group, FC)) +
    geom_col(data=tmp.er, aes(x=group, y=ave_FC, fill=group), width = 0.8) +
    geom_errorbar(data=filter(tmp.er,group!="Monocyte"), aes(
      x=group, y=ave_FC, 
      ymin=ave_FC-se, ymax=ave_FC+se
    ), width=.5, position=position_dodge(.9)) +
    stat_compare_means(label = "p.signif", method = "t.test", size = 3, label.x = 1.8, label.y.npc = rep(0.9,nrow(tmp.er))) +
    scale_y_continuous(limits = clim, sec.axis = sec_axis(~., breaks = c(0,1.0))) +
    ylab("Relative RNA Expression") + xlab("Target Gene") +
    facet_wrap(vars(gene), nrow = 3, scales = "free_y") +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    ggtitle("Monocyte") +
    theme(axis.text.x=element_blank(), 
          axis.text=element_text(size=6), 
          strip.text=element_text(size=6),
          legend.title = element_blank()
    )
  ggsave(width = 8, height = 4, units = "in", device = "png", name, dpi = 300, path = "./results/group3/")
}

# metaplot

allplot <- function(choice)
{
  if(0 %in% choice)
    {for(type in c("TAM","M0"))
    {for(gene in marker$gene)
    {plot1sg(type,gene)}}}
  if(1 %in% choice)
    {for(gene in marker$gene)
    {plot1pl(gene)}}
  if(2 %in% choice)
    {for(con in cond1$condition)
    {plot2(con)}}
}
