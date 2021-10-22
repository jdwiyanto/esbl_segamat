#######################
# ESBL WGS manuscript
#######################
library(tidyverse)
library(janitor)

path = 'G:/My Drive/manuscript1_wgs_esbl/'
setwd(path)

# import community metadata
df = read.csv('metadata_community_raw.csv') # raw metadata

df$age = df$age %>% as.numeric()
df$gender = df$gender %>% str_to_lower %>% factor() 
df$ethnicity = df$ethnicity %>% factor() %>% recode(Native = 'Jakun')
df$education = df$education %>% tolower() %>% factor(levels = c('no', 'dnc_primary', 'primary', 'pmr', 'spm', 'diploma', 'degree'))

chisq.test(table(df$gender))
################
# DEMOGRAPHICS #
################
df$agecat = df$age %>% gtools::quantcut(q=4)
df$bmicat = df %>% mbiome::jd_bmi_cat('bmi') 

demovar = c('agecat', 'gender', 'bmicat', 'ethnicity', 'occupation', 'education')

# convert occupation minor job to others
df2= df %>% mutate(occupation = ifelse(.$occupation == 'agricultural', 'agricultural', ifelse(
  .$occupation == 'children', 'children', ifelse(
    .$occupation == 'homemaker', 'homemaker', ifelse(
      .$occupation =='service', 'service', ifelse(
        .$occupation =='unemployed', 'unemployed', 'others')))))) 

table1 = lapply(demovar, function(x) {
  df2 %>% filter(esbl==T) %>% pull(x) %>% tabyl(show_na=F) %>% data.frame %>% rename('values'='.') %>% add_column(factor=x, .before = 1)
}) %>% bind_rows()

table2 = lapply(demovar, function(x) {
  df2 %>% filter(esbl==F) %>% pull(x) %>% tabyl(show_na=F) %>% data.frame %>% rename('values'='.') %>% add_column(factor=x, .before = 1)
}) %>% bind_rows() %>% select(!factor)

table3 = inner_join(table1, table2, by='values', suffix=c('_pos', '_neg'))

# risk factor analysis
library(lme4)
library(boot)

df.logit = df %>% select(hh, demovar, esbl) %>% na.omit()

# convert occupation minor job to others
df.logit = df.logit %>% mutate(occupation = ifelse(.$occupation == 'agricultural', 'agricultural', ifelse(
  .$occupation == 'children', 'children', ifelse(
    .$occupation == 'homemaker', 'homemaker', ifelse(
      .$occupation =='service', 'service', ifelse(
        .$occupation =='unemployed', 'unemployed', 'others')))))) 

logit.list = lapply(demovar, function(x) {
  logit1 = glmer(esbl ~ get(x) + (1|hh), data = df.logit, family=binomial(link='logit'))
  logit0 = glmer(esbl ~ 1 + (1|hh), data = df.logit, family=binomial(link='logit'))
  
  anova(logit0, logit1)
  summary(logit1)
  
  test = inv.logit(fixef(logit1)) %>% data.frame %>% slice(-1) %>% rownames_to_column() # generate %
  test2 = inv.logit(confint.merMod(logit1, oldNames = F)) %>% data.frame %>% slice(-(1:2)) %>% rownames_to_column() # 95% CI
  test3 = inner_join(test, test2, by='rowname')
})

names(logit.list) = demovar
logit.list2 = logit.list %>% bind_rows(.id = 'factor')
logit.list2$rowname = gsub("get\\(x\\)", '', logit.list2$rowname)

logit.list2


logit.pval = sapply(demovar, function(x) {
  logit1 = glmer(esbl ~ get(x) + (1|hh), data = df.logit, family=binomial(link='logit'))
  logit0 = glmer(esbl ~ 1 + (1|hh), data = df.logit, family=binomial(link='logit'))
  
  anova(logit0, logit1) %>% data.frame %>% slice(-1) %>% pull(Pr..Chisq.)
  
})

logit.pval2 = logit.pval %>% data.frame

table3 # demographics table
table3$values = gsub(',', '-', table3$values) # convert comma to hyphen to avoid csv confusion
write.csv(table3, file = 'table1_demographics.csv', quote=F, row.names = F)

logit.list2 # rate and 95%CI table
logit.list2$rowname = gsub(',', '-', logit.list2$rowname) # convert comma to hyphen to avoid csv confusion
write.csv(logit.list2, file = 'table1_demographics_95CI.csv', quote=F, row.names = F)

logit.pval2 # LRT test
write.csv(logit.pval2, file='table1_demographics_lrt.csv', quote=F, row.names = T)

#########################################################################
# ISOLATION AND PHENOTYPIC PROFILING OF ESBL-PRODUCING ESCHERICHIA COLI #
#########################################################################
library(tidyverse)
library(lme4)
library(boot)
library(viridis)

# import metadata
df.raw = read.csv('metadata_community_raw.csv', header = T)
df.raw$esbl = df.raw$esbl %>% tolower %>% factor

# calculate ESBL carriage rate
model.esbl = glmer(esbl ~ 1 + (1|hh), data = df.raw, family = binomial(link = 'logit'))
summary(model.esbl)
inv.logit(fixef(model.esbl))
inv.logit(confint.merMod(model.esbl, oldNames = F))

# heatmap of AST
df.ast = read.csv('metadata_community_ast.csv', header = T) %>% mutate(across(.cols = matches('[1-9]'), ~recode(., 'I' = 'ns', 'R' = 'ns', 'S' = 's')))
df.ast$esbl = paste0('ESBL-', df.ast$esbl %>% str_to_sentence())

df.ast %>% pivot_longer(!c(household, stoolcode, esbl)) %>% 
  ggplot(aes(x = stoolcode %>% as.factor(), y = name, fill = value)) + 
  geom_tile() + 
  facet_wrap(~esbl, nrow = 1, scales = 'free') +
  ggsci::scale_fill_d3(name = 'Susceptibility', labels = c('Not-susceptible', 'Susceptible')) +
  xlab('Individuals') +
  ylab('Antibiotics') +
  theme_classic() + theme(axis.text.x = element_blank(), legend.position = 'top') 
  
ggsave('fig.ast.pdf', units='in', height = 8, width = 8)  
ggsave('fig1.tiff', units='in', height=7, width=7, compression = 'lzw', dpi=300)

# investigate MDR status
library(tidyverse)
df.mdr = df.ast %>% mutate(across(.cols = !c(household, stoolcode, esbl), ~ifelse(. == 'ns', 1, 0))) %>% 
  mutate(class_aminoglycoside = ifelse(AK30 == 1, 1, 0),
         class_penicillin_inh = ifelse(SAM20 == 1 | TZP110 == 1, 1, 0),
         class_cephalosporin_1 = ifelse(KZ30 == 1, 1, 0),
         class_cephalosporin_4 = ifelse(FEP30 == 1 | CTX30 == 1 | CAZ30 == 1, 1, 0),
         class_quinolone = ifelse(CIP5 == 1 | NA30 == 1, 1, 0),
         class_carbapenem = ifelse(IPM10 == 1, 1 , 0),
         class_nitrofurantoin = ifelse(F300 == 1,1,0),
         class_sulfonamide = ifelse(SXT25==1,1,0),
         class_tetracycline = ifelse(TE30 == 1,1,0)) 

df.mdr %>% select(contains('class')) %>% rowSums %>% table

df.mdr$sum = df.mdr %>% select(contains('class')) %>% rowSums 
mdr.esblpos = df.mdr %>% filter(esbl == 'ESBL-Positive') %>% pull(sum)
mdr.esblneg = df.mdr %>% filter(esbl == 'ESBL-Negative') %>% pull(sum)
mdr.esblpos %>% sd
mdr.esblneg %>% sd
t.test(mdr.esblpos, mdr.esblneg)
##############################################################
# WHOLE-GENOME SEQUENCING OF ESBL-PRODUCING ESCHERICHIA COLI #
##############################################################
library(janitor)
library(viridis)
# srst2 output
srst = read.table('srst2_compiled_rev__compiledResults.txt', sep = '\t', header = T, row.names = 1)

# names of clinical isolate
name.clin = paste0('JD-', c(14:18, 38:40))

# srst df for community only
srst.com = srst %>% rownames_to_column() %>% filter(!(rowname %in% name.clin)) %>% column_to_rownames()

# MLST
library(viridis)
srst.com$ST %>% tabyl %>% 
  dplyr::arrange(desc(n)) %>% 
  ggplot(aes(x = reorder(., n), y = n, fill = .)) + 
  geom_col() + 
  coord_flip() + 
  theme_bw() +
  scale_fill_viridis(discrete = T) + 
  theme(legend.position = 'none') +
  xlab('Strain Type') +
  ylab('No. isolate')

ggsave('fig.mlst.tiff', units='in', height = 8, width = 8, dpi=300)  

# no of unique st
srst.com$ST %>% unique %>% length

# check for intrafamilial transmission
srst.com %>% dplyr::select(ST) %>% filter(ST == 131)

srst.arg = read.table('srst2_argannot_compiled__compiledResults.txt', header = T, sep = '\t') # ARGannot3 database
srst.vf = read.table('srst2_VF_compiled__compiledResults.txt', header = T, sep = '\t') # VF database
srst.pf = read.table('srst2_plasmidfinder_compiled__compiledResults.txt', header = T, sep = '\t') # PlasmidFinder database
srst.pr = read.table('srst2_plasmidreplicon_compiled__compiledResults.txt', header = T, sep = '\t') # Plasmid18Replicons database

wgs.code = read.csv('metadata_decoder.csv', header = T) %>% filter(note == 'community') %>% pull(stoolcode) # pull stoolcodes of WGS isolate 
wgs.newcode = read.csv('metadata_decoder.csv', header = T) %>% filter(note == 'community') %>% arrange((as.numeric(stoolcode))) %>% pull(new_code) # pull new code for WGS isolate

srst.all = bind_rows(
            srst.arg %>% pivot_longer(!Sample) %>% mutate(value2 = ifelse(value == '-', 'absent', 'present'), db = 'Antibiotic Resistance Gene'), 
            srst.vf %>% pivot_longer(!Sample) %>% mutate(value2 = ifelse(value == '-', 'absent', 'present'), db = 'Virulence Finder'), 
            srst.pf %>% pivot_longer(!Sample) %>% mutate(value2 = ifelse(value == '-', 'absent', 'present'), db = 'Plasmid Finder') 
            # srst.pr %>% pivot_longer(!Sample) %>% mutate(value2 = ifelse(value == '-', 'absent', 'present'), db = 'Plasmid Replicon')
            )

# add setting information
srst.all = srst.all %>% mutate(setting = ifelse(Sample %in% name.clin, 'clinical', 'community'))

# save srst.all as supplementary information to PLOS one
srst.all.rev = srst.all
colnames(srst.all.rev) = c('Sample', 'Gene', 'Variant', 'Presence', 'Database', 'Setting')
write.csv(srst.all.rev, file='draft4_submission/revision1/resistome_list.csv', quote=F, row.names = F)


# virulence factor correlation with setting?
srst.vf.com2
srst.vf2 = srst.vf %>% column_to_rownames('Sample') %>% mutate(across(.fns = ~ifelse(.=='-', 0, 1))) 

srst.vf3 = srst.vf2 %>% rownames_to_column('sample') %>% 
  mutate(setting = ifelse(sample %in% name.clin, 'clinical', 'community')) %>% 
  column_to_rownames('sample') 

vf.com = (srst.vf3 %>% filter(setting == 'community') %>% select(!setting) %>% colSums) / 32 * 100
vf.clin = (srst.vf3 %>% filter(setting == 'clinical') %>% select(!setting) %>% colSums) / 8 * 100

bind_rows(
data.frame(vf.com, vf.clin) %>% rownames_to_column('sample') %>% 
  mutate(value2 = vf.com - vf.clin) %>% 
  arrange(-value2) %>%
  slice_head(n=20)
,
data.frame(vf.com, vf.clin) %>% rownames_to_column('sample') %>% 
  mutate(value2 = vf.com - vf.clin) %>% 
  arrange(-value2) %>%
  slice_tail(n=20)
) %>% mutate(dir = ifelse(value2 > 0, 'community', 'clinical'),
             sample = paste0('*', sample, '*')) %>% 
  ggplot(aes(x=reorder(sample, value2), y=value2, fill=dir)) + geom_col() + coord_flip() +
  theme_bw() +
  labs(y = 'Proportion difference (%)',
       x = 'Virulence Gene') +
  scale_fill_manual(name = '', 
                    values = c('steelblue', 'orange'),
                    labels = c('Higher in clinical isolates',
                               'Higher in community isolates')) +
  theme(legend.position = 'top',
        axis.text.y = element_markdown())

ggsave('draft4_submission/revision1/fig_vf_diff.pdf', height=6, width=8)

# how many VF are exclusiver to one setting?
vf.unq.com = srst.vf3 %>% filter(setting == 'community') %>% rownames_to_column() %>% select(!setting) %>% pivot_longer(!rowname) %>% filter(value==1) %>% pull(name) %>% unique
vf.unq.clin = srst.vf3 %>% filter(setting == 'clinical') %>% rownames_to_column() %>% select(!setting) %>% pivot_longer(!rowname) %>% filter(value==1) %>% pull(name) %>% unique 

(vf.unq.com %in% vf.unq.clin)%>% table
(vf.unq.clin %in% vf.unq.com)%>% table
bind_rows(
vf.unq.com %>% data.frame %>% mutate(setting = 'community'),
vf.unq.clin %>% data.frame %>% mutate(setting = 'clinical')) %>% arrange(`.`)

rowSums(srst.vf2)
vf.ttest = srst.vf2 %>% rownames_to_column('Sample') %>% 
  mutate(sum = vf.sum,
         setting = ifelse(Sample %in% name.clin, 'clinical', 'community')) %>% select(sum, setting) 

t.test()

################################
# basic description of gene data
################################

# arg database
df.count.arg = srst.arg %>% 
  mutate(across(.cols = !Sample, ~ifelse(.=="-", 0, 1))) %>% 
  mutate(setting = ifelse(Sample %in% name.clin, 'clinical', 'community')) 

# plasmid database
df.count.pf = srst.pf %>% mutate(across(.cols = !Sample, ~ifelse(.=="-", 0, 1))) %>% 
  mutate(setting = ifelse(Sample %in% name.clin, 'clinical', 'community')) 

# virulence factor database
df.count.vf = srst.vf %>% mutate(across(.cols = !Sample, ~ifelse(.=="-", 0, 1))) %>% 
  mutate(setting = ifelse(Sample %in% name.clin, 'clinical', 'community')) 

df.count = df.count.arg

# no of unique genes
df.count %>% select(!c(Sample, setting)) %>% ncol

# no unique gene classes
df.count %>% select(!c(Sample, setting)) %>% colnames() %>% str_split_fixed(pattern = '_', n = 2) %>% data.frame %>% select(X2) %>% unique() %>% nrow # gene class

# gene frequency
apply(df.count %>% select(!c(Sample, setting)), 2, sum) %>% sort() # gene frequency

# gene average
df.count %>% select(!c(Sample, setting)) %>% rowSums %>% summary() # avg gene carriage
df.count %>% select(!c(Sample, setting)) %>% rowSums %>% sd()

# any difference between mean carriage of community v clinical isolate?
# use three databases: df.count.arg, df.count.vf, df.count.pf
t.clin = df.count.vf %>% filter(Sample %in% name.clin) %>% select(!c(Sample, setting)) %>% rowSums
t.comm = df.count.vf %>% filter(!(Sample %in% name.clin)) %>% select(!c(Sample, setting)) %>% rowSums
t.test(t.clin, t.comm)

##################################################
# differences in gene carriage across ctxm variant
##################################################

lapply(c("CTX-M-15", "CTX-M-55", "CTX-M-3", "CTX-M-69", "CTX-M-27", "CTX-M-65"), function(x) {
  df.count.vf %>% select(!setting) %>% inner_join(df.ctxm, by='Sample') %>% 
  filter(CTXM == x) %>%
  select(!c(Sample, CTXM)) %>%
  rowSums %>% summary
})
 
lapply(c("CTX-M-15", "CTX-M-55", "CTX-M-3", "CTX-M-69", "CTX-M-27", "CTX-M-65"), function(x) {
  df.count.arg %>% select(!setting) %>% inner_join(df.ctxm, by='Sample') %>% 
    filter(CTXM == x) %>%
    select(!c(Sample, CTXM)) %>%
    rowSums %>% sd
})

 



#logit df for ARF v antibiotic susceptibility test result

# AST profile for community isolates
df.ast.wgs = df.ast %>% filter(stoolcode %in% wgs.code) %>% mutate(newcode = paste0('JD-', wgs.newcode)) %>% dplyr::select(!c(household, stoolcode, esbl)) %>%
  mutate_at(vars(1:13), ~ifelse(. == 's', 0, 1)) # filter AST result to only WGS isolate

# Antibiotic resistance gene profile for community isolates
srst.arg.com = srst.arg %>% filter(Sample %in% df.ast.wgs$newcode) # filter srst.arg result to only WGS community isolate
srst.arg.com = srst.arg.com[order(match(srst.arg.com$Sample, df.ast.wgs$newcode)),] %>% 
  remove_rownames() %>% column_to_rownames(var = 'Sample') %>% 
  mutate_at(vars(everything()), ~ifelse(.=='-', 0, 1)) # reorder srst.arg.com to follow ast df order

# plasmid carriage profile for community isolates
srst.pf.com

#############################################################
# logit for antibiotic resistance genes and plasmid carriage

srst.arg.com = srst.arg %>% filter(Sample %in% srst.pf.com$newcode) %>% mutate(across(.cols = !Sample, ~ifelse(.=='-', 0, 1))) %>% rename(newcode = Sample)
srst.pf.com # plasmid profile

# combine both
arg.pf.df = inner_join(srst.arg.com, srst.pf.com, by='newcode')

# manually check variables with <2 levels
coltokeep = sapply(1:ncol(arg.pf.df), function(x) arg.pf.df[[x]] %>% factor %>% nlevels) != 1

arg.pf.df2 = arg.pf.df[,coltokeep] 
arg.pf.df2 %>% colnames

arg.pf.model = lapply(colnames(arg.pf.df2)[2:31], function(var.ast) {
  lapply(colnames(arg.pf.df2)[32:58], function(var.pf) {
    
    test3 = glm(var.ast %>% get ~ var.pf %>% get, data = arg.pf.df2, family = binomial(link='logit')) %>% summary
    # test3$coefficients
    test4 = test3$coefficients[2,] %>% t %>% data.frame %>% mutate(arg = var.ast, plasmid = var.pf)
    return(test4)
  })
  
}) %>% bind_rows() %>% filter(Pr...z.. < 0.05)

arg.pf.model %>% filter(plasmid == 'I1_Alpha')

write.csv(arg.pf.model, file = 'table.arg.pf.csv', quote=F)

############################################################

##############################################################################
# logit for antibiotic resistance gene versus antibiotic susceptibility profile

srst.arg.com
df.ast.wgs 

# combine both
arg.ast.df = inner_join(srst.arg.com, df.ast.wgs, by='newcode')

# manually check variables with <2 levels
coltokeep = sapply(1:ncol(arg.ast.df), function(x) arg.ast.df[[x]] %>% factor %>% nlevels) != 1

arg.ast.df2 = arg.ast.df[,coltokeep] 
arg.ast.df2 %>% colnames

arg.ast.model = lapply(colnames(arg.ast.df2)[2:31], function(var.ast) {
  lapply(colnames(arg.ast.df2)[32:40], function(var.pf) {
    
    test3 = glm(var.ast %>% get ~ var.pf %>% get, data = arg.ast.df2, family = binomial(link='logit')) %>% summary
    # test3$coefficients
    test4 = test3$coefficients[2,] %>% t %>% data.frame %>% mutate(arg = var.ast, ast = var.pf)
    return(test4)
  })
  
}) %>% bind_rows() %>% filter(Pr...z.. < 0.05)

arg.ast.model %>% count(ast)

write.csv(arg.ast.model, file = 'table.arg.ast.csv', quote=F)

##############################################################################
##########################################################
# logit for susceptibility profile versus plasmid carriage
df.ast.wgs # AST profile
srst.pf.com = srst.pf %>% filter(Sample %in% df.ast.wgs$newcode) %>% mutate(across(.cols = !Sample, ~ifelse(.=='-', 0, 1))) %>% dplyr::rename(newcode = Sample) 
ast.pf.df = inner_join(df.ast.wgs, srst.pf.com, by='newcode')

# manually check variables with <2 levels
coltokeep = sapply(1:ncol(ast.pf.df), function(x) ast.pf.df[[x]] %>% factor %>% nlevels) != 1

ast.pf.df2 = ast.pf.df[,coltokeep] 
ast.pf.df2 %>% colnames

ast.pf.model = lapply(colnames(ast.pf.df2)[1:9], function(var.ast) {
  lapply(colnames(ast.pf.df2)[11:37], function(var.pf) {
    
    test3 = glm(var.ast %>% get ~ var.pf %>% get, data = ast.pf.df2, family = binomial(link='logit')) %>% summary
    # test3$coefficients
    test4 = test3$coefficients[2,] %>% t %>% data.frame %>% mutate(ast = var.ast, plasmid = var.pf)
    return(test4)
  })
  
}) %>% bind_rows() %>% filter(Pr...z.. < 0.1)

write.csv(ast.pf.model, file = 'table.ast.pf.csv', quote=F)
 
#####################################################################################
# transmission of ESBL-EC between community dwellers and clinical patients in segamat
#####################################################################################
library(tidyverse)

#generate metadata for phandago viewing

df.phylo = bind_rows(
read.csv('metadata_wgs_community.csv', header = T) %>% mutate(id = paste0('JD-', new_code),
                                                              sex = (.$gender %>% str_sub(1,1)),
                                                              source = 'faeces',
                                                              setting = 'community') %>% dplyr::select(c(id, age, sex, setting, source)),
read.csv('metadata_wgs_clinical.csv', header = T) %>% mutate(id = paste0('JD-', id),
                                                             setting = 'clinical') %>% rename(source = isolation_source) %>% 
  dplyr::select(c(id, age, sex, setting, source))
)

# add mlst data
phylo.mlst = srst %>% rownames_to_column(var = 'id') %>% dplyr::select(c(id, ST))
phylo.mlst = phylo.mlst[order(match(phylo.mlst$id, df.phylo$id)),]

df.phylo$st = phylo.mlst$ST
write.csv(df.phylo, 'metadata_phandango.csv', quote = F, row.names = F)

# phylogenetic tree
library(ggtree)
library(ggnewscale)
library(ggsci)

phylo.tree = read.tree(file = 'my_tree.newick')

# change tip label name based on setting
d1 = data.frame(id = phylo.tree$tip.label,
                Setting = df.phylo[order(match(df.phylo$id, phylo.tree$tip.label)),]$setting %>% factor) %>% 
  mutate(across(.cols='Setting', str_to_sentence)) %>% 
  column_to_rownames(var = 'id')

p = ggtree(phylo.tree) + 
  geom_tiplab(align=T, linesize=.7) + 
  geom_hilight(node=62, fill='cyan')
  #geom_cladelabel(node=62, label='ST131', offset=0.003, color='red2', align=T)

library(ggsci); library(viridis)
p1 = gheatmap(p, d1, offset = 0.006, width = 0.2, font.size = 2, hjust = 0.5) + scale_fill_aaas(name='Setting')

d2 = data.frame(id = phylo.tree$tip.label,
                MLST = phylo.mlst[order(match(phylo.mlst$id, phylo.tree$tip.label)),]$ST) %>% column_to_rownames(var = 'id')

library(ggnewscale)
p2 = p1 + new_scale_fill()
p3 = gheatmap(p2, d2, offset=.013, width=0.2, font.size=2, hjust=0.5) + scale_fill_igv(name='MLST')

# add ctx-m annotation
df.arg.tree = srst.arg[order(match(srst.arg$Sample, rownames(d1))),] %>% select(Sample, contains('CTX')) 
df.arg.tree$CTX.M.1_Bla = df.arg.tree$CTX.M.1_Bla %>% str_split_fixed(pattern='_', n=2) %>% data.frame %>% pull(X1)
df.arg.tree$CTX.M.9_Bla = df.arg.tree$CTX.M.9_Bla %>% str_split_fixed(pattern='_', n=2) %>% data.frame %>% pull(X1)
colnames(df.arg.tree) = colnames(df.arg.tree) %>% str_split_fixed(pattern='_', n=2) %>% data.frame %>% pull(X1)

df.arg.tree2 = df.arg.tree %>% mutate("CTX-M-15" = ifelse(CTX.M.1 == 'CTX-M-15', 'Present', 'Absent'),
                       "CTX-M-55" = ifelse(CTX.M.1 == 'CTX-M-55', 'Present', 'Absent'),
                       "CTX-M-3" = ifelse(CTX.M.1 == 'CTX-M-3', 'Present', 'Absent'),
                       "CTX-M-69" = ifelse(CTX.M.1 == 'CTX-M-69', 'Present', 'Absent'),
                       "CTX-M-27" = ifelse(CTX.M.9 == 'CTX-M-27', 'Present', 'Absent'),
                       "CTX-M-65" = ifelse(CTX.M.9 == 'CTX-M-65', 'Present', 'Absent')) %>%
  select(-c(CTX.M.1, CTX.M.9)) %>% column_to_rownames(var='Sample')
rownames(df.arg.tree2) = phylo.tree$tip.label

p4 = p3 + new_scale_fill()
p5 = gheatmap(p4, df.arg.tree2, offset=0.020,width=1.2, font.size=1.75, hjust=0.5) + scale_fill_d3(name = 'CTX-M')
p5
ggsave('fig.tree2.pdf', height = 10, width = 19, units='in')
ggsave('Fig4.tiff', units = 'in', height=8.75, width=7.5, dpi=300, compression='lzw')

# pangenomic comparison with publicly available ST131 data
library(tidyverse)
library(ggtree)
library(viridis)
library(ggsci)

df.pub = read.csv('roary_out/metadata.curated.csv', header=T) %>% dplyr::select(run, country, status) %>% mutate(across(c(country, status), str_to_title))
d1.pub = df.pub %>% column_to_rownames(var = 'run') %>% dplyr::select(status) %>% rename(setting = status)
d2.pub = df.pub %>% column_to_rownames(var = 'run') %>% dplyr::select(country) 
tree.pub = read.tree(file = 'roary_out/my_tree.newick')
p.pub = ggtree(tree.pub, branch.length='none', layout='circular') %<+% df.pub + geom_tiplab(aes(color=country), align=T, linetype=NA, size=1.5, offset=8, hjust=0.5)  + geom_tippoint(aes(color=country))+ scale_color_igv(name='Country')

p.pub.labeled = p.pub + geom_hilight(node=369, fill='cyan') + geom_highlight(node=286, fill='magenta')
p2.pub = gheatmap(p.pub.labeled, d1.pub, colnames_position = 'top', offset=15, hjust=0.5, font.size = 1.5, width = 0.05, colnames_angle = 90, colnames_offset_y = 1) + scale_fill_viridis(discrete=T, name='Setting')
p3.pub = p2.pub + new_scale_fill()
gheatmap(p3.pub, d2.pub, colnames_position = 'top', offset=18, hjust=0.5, font.size = 1.5, width = 0.05, colnames_angle = 90, colnames_offset_y = 1) + scale_fill_igv(name='Country')

ggsave('fig_st131_public.pdf', units='in', height= 12, width=12)
ggsave('Fig3.tiff', units = 'in', height = 7.5, width=7.5, dpi=300, compression='lzw')

# CTXM variants
library(ggtext)
library(ggsci)

fig.com.ctxm = srst.com %>% dplyr::select(contains('CTX.M')) %>% rownames_to_column() %>% pivot_longer(!rowname) %>% dplyr::filter(value != '-') %>% count(value, name) %>% mutate(xtitle = paste0(value, '<br>(n=',n, ')')) %>% ggplot(aes(x=reorder(xtitle, n), y=n, fill=name)) + geom_col() + ylim(c(0,10)) + coord_flip() + facet_wrap(~name, scales='free')  + scale_fill_d3() + theme_bw()+ theme(legend.position = 'none', axis.text.y = element_markdown()) + xlab('CTX-M Variants') + ylab('Count')

ggsave(filename = 'fig.com.ctxm.tiff', units = 'in', plot = fig.com.ctxm, height = 8, width=8, dpi=300)


srst.com %>% dplyr::select(contains('CTX.M')) %>% rownames_to_column() %>% pivot_longer(!rowname) %>% dplyr::filter(value != '-') %>% count(value, name) 

# revealed ctxm-9 variant differences between cluster1 and cluster2
srst %>% select(contains('CTX')) %>% rownames_to_column() %>% filter(rowname %in% cl1.name)
srst %>% select(contains('CTX')) %>% rownames_to_column() %>% filter(rowname %in% cl2.name)

srst %>% select(contains('CTX')) %>% rownames_to_column() %>% filter(rowname %in% name.clin)
####################################################################################################################
# master figure annotating across MLST, phenotypic resistance, antibiotic resistance gene carried, plasmid group carried
####################################################################################################################

srst.st.com = srst.com %>% select(ST) %>% rownames_to_column(var='newcode') # community mlst data
df.ast.wgs # community Ast data
srst.arg.com # community arg data
srst.pf.com # community plasmid data

cordf = inner_join(df.ast.wgs, srst.arg.com,by='newcode') %>% inner_join(srst.pf.com, by='newcode') %>% column_to_rownames(var='newcode') %>% select(!c(KZ30, IPM10, TZP110, CTX30, AMPH_Ecoli_Bla, AmpC2_Ecoli_Bla, ErmC_MLS, MrdA_Bla, Q1)) 
colnames(cordf)[1:9] = paste0('AST_', colnames(cordf[1:9]))
colnames(cordf)[10:39] = paste0('ARG_', colnames(cordf[10:39]))
colnames(cordf)[40:66] = paste0('PF_', colnames(cordf[40:66]))


cordf.cor = cordf %>% cor
cordf.pval = cor.mtest(cordf, conf.level = 0.95)

library(corrplot)

tiff(filename='fig_corrplot.tiff', height=9, width=9, units='in', res=300)
corrplot(cordf.cor, 
         p.mat = cordf.pval$p, 
         sig.level = 0.05,
         insig = 'blank',
         pch.col = 'grey20',
         pch.cex=1,
         method='circle', 
         order='hclust',
         type='lower', 
         tl.cex = 0.7, 
         tl.col = 'dodgerblue4', 
         tl.srt=45)
dev.off()

# turn into data.frame

cordf.cor
rownames(cordf.pval$p) = rownames(cordf.cor)
colnames(cordf.pval$p) = rownames(cordf.cor)

rownames(cordf.pval$lowCI) = rownames(cordf.cor)
colnames(cordf.pval$lowCI) = rownames(cordf.cor)

rownames(cordf.pval$uppCI) = rownames(cordf.cor)
colnames(cordf.pval$uppCI) = rownames(cordf.cor)

cordf.cor[lower.tri(cordf.cor)] = NA
cordf.pval$p[lower.tri(cordf.pval$p)] = NA
cordf.pval$lowCI[lower.tri(cordf.pval$lowCI)] = NA
cordf.pval$uppCI[lower.tri(cordf.pval$uppCI)] = NA

cor.cor = na.omit(reshape2::melt(cordf.cor)) %>% rename(cor=value) %>% mutate(pair = paste0(Var1, '_', Var2)) %>% select(pair, Var1, Var2, cor)
cor.p = na.omit(reshape2::melt(cordf.pval$p)) %>% rename(pval = value) %>% mutate(pair = paste0(Var1, '_', Var2)) %>% select(pair, pval)
cor.low = na.omit(reshape2::melt(cordf.pval$lowCI)) %>% rename(lowci = value) %>% mutate(pair = paste0(Var1, '_', Var2)) %>% select(pair, lowci)
cor.hi = na.omit(reshape2::melt(cordf.pval$uppCI)) %>% rename(uppci = value) %>% mutate(pair = paste0(Var1, '_', Var2)) %>% select(pair, uppci)

cor.df.final = inner_join(cor.cor, cor.p, by='pair') %>% inner_join(cor.low, by='pair') %>% inner_join(cor.hi, by='pair') %>% mutate(padj = p.adjust(pval, method = 'BH')) %>% 
  mutate(dupli = ifelse(Var1 == Var2, 'dupli', 'unique')) %>% filter(dupli=='unique') %>% 
  mutate(pair1 = str_split_fixed(Var1, pattern='_', n=2) %>% data.frame %>% pull(X1),
         pair2 = str_split_fixed(Var2, pattern='_', n=2) %>% data.frame %>% pull(X1),
         pair3 = paste0(pair1, '_', pair2)) %>%
  filter(padj < 0.05) %>% mutate(dir = ifelse(cor<0, 'neg', 'pos')) 
  
cor.df.final %>% ggplot(aes(x=pair, y=cor, fill=dir)) + geom_col() + coord_flip() + theme_bw() + facet_grid(~pair3, scales='free_y') + ylim(c(-1,1)) + 
  scale_fill_manual(values=c('indianred3', 'dodgerblue2')) + geom_hline(yintercept = 0, linetype='dashed') + theme(legend.position = 'none') +
  ylab('Correlation coefficient') + xlab('Pair') + geom_errorbar(aes(ymin = lowci, ymax = uppci, width=0.1))

ggsave('fig_corr.tiff', height=16,width=8, dpi=300, units = 'in')  

cor.df.final %>% arrange(-abs(cor))
cor.df.final %>% count(pair3)

cor.df.final %>% filter(pair3 == 'ARG_PF') %>% pull(Var2) %>% tabyl
save.image()

##########################
# complex heatmap analysis
##########################
library(ComplexHeatmap)

cordf2 = inner_join(df.ast.wgs, srst.arg.com,by='newcode') %>% 
  inner_join(srst.pf.com, by='newcode') %>% 
  column_to_rownames(var='newcode') 

rownames(cordf2) = gsub('JD', 'community', rownames(cordf2))
colnames(cordf2) = colnames(cordf2) %>% str_split_fixed(pattern='_', n=2) %>% data.frame %>% pull(X1)

heatmap.col = structure(c(4,7), names=c('0', '1'))

# tiff('fig.heatmap.vf.tiff', height=8,width=24, res=300, units = 'in')
Heatmap(cordf2, 
        name = 'Legend',
        column_title_side = 'top',
        col = heatmap.col,
        column_names_rot = 90, 
        column_names_gp = gpar(fontsize=8),
        row_title = 'Isolate', 
        row_title_side = 'right', 
        row_title_rot = 90, 
        row_km=2,
        row_gap = unit(0.5, 'cm'), 
        row_names_gp = gpar(fontsize=8),
        cluster_rows = T, 
        row_dend_width = unit(3, 'cm'),
        cluster_columns = T, 
        show_column_dend = F, 
        column_split = factor(c(rep('AST', 13), 
                                rep('ARGannot3', 34), 
                                rep('PlasmidFinder', 28))), 
        column_gap = unit(1, 'cm'),
        heatmap_legend_param = list(title='', at = c(0,1),
                                    labels = c('Negative', 'Positive')
        ))
# dev.off()


save.image()

# no plasmid and arg carried by each CTXM

cordf3 %>% mutate(stat = ifelse(CTX.M.1_Bla == 1, 'CTXM1', 'CTXM9')) %>% select(1:28, stat) %>% pivot_longer(!stat) %>% filter(value!=0) %>% group_by(stat) %>% summarise(n = n()) 

cordf3 %>% mutate(stat = ifelse(CTX.M.1_Bla == 1, 'CTXM1', 'CTXM9')) %>% select(29:62, stat) %>% pivot_longer(!stat) %>% filter(value!=0) %>% group_by(stat) %>% summarise(n = n()) 

# identify cluster link with geography

cluster1 = c(9,20,24,27,31,32,7,28,8,36,30,22,25,21,23,1,34,33)
cluster2 = c(6,29,26,19,12,2,13,37,5,10,11,35,4,3)

cl1.name = paste0('JD-', cluster1)
cl2.name = paste0('JD-', cluster2)
cl.all.name = c(cl1.name, cl2.name)

df.ctxm = df.arg.tree %>% mutate(CTXM = ifelse(CTX.M.1 != '-', CTX.M.1, CTX.M.9)) %>% select(Sample, CTXM)

library(tidyverse)
df.decoder = read.csv('metadata_decoder.csv', header=T)
df.decoder2 = df.decoder %>% filter(new_code %in% c(cluster1, cluster2))
newstoolcode = df.decoder %>% filter(new_code %in% c(cluster1, cluster2)) %>% pull(stoolcode)

df.demo = read.csv('C:\\Users\\iajac\\Dropbox\\PhD\\questionnaire\\qs_demography_extended.csv', header=T)
df.demo2 = df.demo %>% filter(stoolcode %in% newstoolcode)
df.demo3 = df.demo2[order(match(df.demo2$stoolcode, df.decoder2$stoolcode)),] %>% mutate(newcode = df.decoder2$new_code)

df.demo4 = df.demo3 %>% mutate(Sample = paste0('JD-', newcode)) %>% inner_join(df.ctxm, by='Sample')

# plot geo distance
segamat.lat = df.demo4 %>% pull(coordinate) %>% str_split_fixed(pattern =',', n=2) %>% data.frame %>% pull(X1)
segamat.lon = df.demo4 %>% pull(coordinate) %>% str_split_fixed(pattern =',', n=2) %>% data.frame %>% pull(X2)

library(ggsci)
df.demo4 %>% mutate(cluster = ifelse(newcode %in% cluster1, 'blaCTX-M-1', 'blaCTX-M-9'),
                    lat = segamat.lat, 
                    lon = segamat.lon) %>% 
  select(lat, lon, cluster, CTXM) %>% mutate(across(!c(CTXM,cluster), as.numeric)) %>% 
  ggplot(aes(x=lon, y=lat, color=CTXM)) + geom_point(size=2.5, show.legend=T) + coord_quickmap() + theme_bw() + labs(x = 'Longitude', y = 'Latitude') + theme(legend.position = 'top') + scale_color_aaas(name = 'CTX-M variant') + facet_wrap(~cluster)

ggsave('fig_geo.pdf', units = 'in', height=8,width=8)


################################################
# complex heatmap with clinical isolate together
################################################

df.clin = read.csv('metadata_wgs_clinical.csv', header=T)
df.clin.ast = df.clin %>% mutate(id = paste0('JD-', id)) %>% select(!c(barcode, age, sex, diagnosis, isolation_date, isolation_source)) %>% rename(Sample=id) %>% mutate(across(!Sample, ~ifelse(.=='S', 0, 1))) %>% select(!c(FOX30, TZP110)) %>% column_to_rownames(var='Sample')

df.ast.wgs2 = df.ast.wgs %>% column_to_rownames('newcode')

mutual.ast = colnames(df.clin.ast ) %in% colnames(df.ast.wgs2)
mutual.ast2 = colnames(df.ast.wgs2 ) %in% colnames(df.clin.ast)

df.clin.ast2 = df.clin.ast[,mutual.ast]
df.ast.wgs3 = df.ast.wgs2[,mutual.ast2]

df.ast.all = bind_rows(df.ast.wgs3, df.clin.ast2) %>% rownames_to_column(var='Sample')
cordf3 # combination of clinical and community segamat isolates 

cordf.all = inner_join(cordf3 %>% rownames_to_column(var='Sample'), df.ast.all, by='Sample') %>% column_to_rownames(var='Sample')

# beta-lactamase
cordf.bla = cordf.all %>% select(contains('Bla')) %>% rownames_to_column(var='Sample')
colnames(cordf.bla) = str_split_fixed(colnames(cordf.bla), pattern='_', n=2) %>% data.frame %>% pull(X1)

# aminoglycoside
cordf.agly =  cordf.all %>% select(contains('Agly')) %>% rownames_to_column(var='Sample')
colnames(cordf.agly) = str_split_fixed(colnames(cordf.agly), pattern='_', n=2) %>% data.frame %>% pull(X1)

# phe
cordf.phe =  cordf.all %>% select(contains('_Phe')) %>% rownames_to_column(var='Sample')
colnames(cordf.phe) = str_split_fixed(colnames(cordf.phe), pattern='_', n=2) %>% data.frame %>% pull(X1)

# tmt
cordf.tmt =  cordf.all %>% select(contains('_tmt')) %>% rownames_to_column(var='Sample')
colnames(cordf.tmt) = str_split_fixed(colnames(cordf.tmt), pattern='_', n=2) %>% data.frame %>% pull(X1)

# mls
cordf.MLS =  cordf.all %>% select(contains('_MLS')) %>% rownames_to_column(var='Sample')
colnames(cordf.MLS) = str_split_fixed(colnames(cordf.MLS), pattern='_', n=2) %>% data.frame %>% pull(X1)

# Fcyn
cordf.Fcyn =  cordf.all %>% select(contains('_Fcyn')) %>% rownames_to_column(var='Sample')
colnames(cordf.Fcyn) = str_split_fixed(colnames(cordf.Fcyn), pattern='_', n=2) %>% data.frame %>% pull(X1)

#colistin
cordf.colistin =  cordf.all %>% select(contains('_colistin')) %>% rownames_to_column(var='Sample')
colnames(cordf.colistin) = str_split_fixed(colnames(cordf.colistin), pattern='_', n=2) %>% data.frame %>% pull(X1)

# flq
cordf.flq =  cordf.all %>% select(contains('_flq')) %>% rownames_to_column(var='Sample')
colnames(cordf.flq) = str_split_fixed(colnames(cordf.flq), pattern='_', n=2) %>% data.frame %>% pull(X1)

# sul
cordf.sul =  cordf.all %>% select(contains('_sul')) %>% rownames_to_column(var='Sample')
colnames(cordf.sul) = str_split_fixed(colnames(cordf.sul), pattern='_', n=2) %>% data.frame %>% pull(X1)

# tet
cordf.tet =  cordf.all %>% select(contains('_tet')) %>% rownames_to_column(var='Sample')
colnames(cordf.tet) = str_split_fixed(colnames(cordf.tet), pattern='_', n=2) %>% data.frame %>% pull(X1)

cordf.ast = cordf.bla %>% inner_join(cordf.agly, by='Sample') %>%
  inner_join(cordf.phe, by='Sample') %>% 
  inner_join(cordf.tmt, by='Sample') %>% 
  inner_join(cordf.MLS, by='Sample') %>% 
  inner_join(cordf.Fcyn, by='Sample') %>% 
  inner_join(cordf.colistin, by='Sample') %>% 
  inner_join(cordf.flq, by='Sample') %>% 
  inner_join(cordf.sul, by='Sample') %>% 
  inner_join(cordf.tet, by='Sample') %>% select(!c(Sample, CTX.M.1, CTX.M.9))

# annotate sample name based on setting
colnames(cordf.all) =colnames(cordf.all) %>% str_split_fixed(pattern='_', n=2) %>% data.frame %>% pull(X1)

cordf.all2 = bind_rows(
  cordf.all %>% filter(setting == 1) %>% 
    rownames_to_column(var='rowname') %>% 
    mutate(rowname = paste0('community-', str_split_fixed(rowname, pattern='-', n=2) %>% data.frame %>% pull(X2))) %>% 
    column_to_rownames(),
  
  cordf.all %>% filter(setting == 0) %>% 
    rownames_to_column(var='rowname') %>% 
    mutate(rowname = paste0('clinical-', str_split_fixed(rowname, pattern='-', n=2) %>% data.frame %>% pull(X2))) %>%
    column_to_rownames()
) %>% select(!setting)

cordf.all3 = cordf.all2 %>% mutate(Sample = paste0("JD-", str_split_fixed(rownames(cordf.all2), pattern = '-', n=2) %>% data.frame %>% pull(X2))) %>% inner_join(df.ctxm, by='Sample') %>% 
  select(!c(Sample, CTX.M.1, CTX.M.9)) %>% 
  mutate(CTX.M.15 = ifelse(CTXM == "CTX-M-15", 1, 0),
         CTX.M.55 = ifelse(CTXM == "CTX-M-55", 1, 0),
         CTX.M.3 = ifelse(CTXM == "CTX-M-3", 1, 0),
         CTX.M.69 = ifelse(CTXM == "CTX-M-69", 1, 0),
         CTX.M.65 = ifelse(CTXM == "CTX-M-65", 1, 0),
         CTX.M.27 = ifelse(CTXM == "CTX-M-27", 1, 0)) %>%
  select(!CTXM)
  
rownames(cordf.all3) = rownames(cordf.all2)

# remove old arg columns from 29:60
cordf.all4 = cordf.all3 %>% select(!c(29:60)) %>% bind_cols(cordf.ast)
colnames(cordf.all4)

# add MLST profile
cordf.mlst = srst %>% rownames_to_column %>% select(rowname, ST) %>% mutate(ST131 = ifelse(ST == '131', 1, 0),
                               ST155 = ifelse(ST == '155', 1, 0),
                               ST69 = ifelse(ST == '69', 1, 0),
                               ST38 = ifelse(ST == '38', 1, 0),
                               ST7401 = ifelse(ST == '7401', 1, 0),
                               ST3580 = ifelse(ST == '3580', 1, 0),
                               ST1421 = ifelse(ST == '1421', 1, 0),
                               ST117 = ifelse(ST == '117', 1, 0),
                               ST361 = ifelse(ST == '361', 1, 0)) %>% select(!ST)

cordf.all5 = cordf.all4 %>% mutate(rowname = paste0("JD-", str_split_fixed(rownames(cordf.all2), pattern = '-', n=2) %>% data.frame %>% pull(X2)))%>% inner_join(cordf.mlst, by='rowname') %>% select(!rowname)
rownames(cordf.all5) = rownames(cordf.all4)

rownames(cordf.all5) = mgsub::mgsub(rownames(cordf.all5), pattern = c('community', 'clinical'), replacement = c('JD', 'JD'))

pdf('draft4_submission/revision1/fig.heatmap.all_rev.pdf', height=9, width=18)
Heatmap(cordf.all5, 
        name = 'Legend',
        column_title_side = 'top',
        col = heatmap.col,
        column_names_rot = 45, 
        column_names_gp = gpar(fontsize=7),
        row_title = 'Isolate', 
        row_title_side = 'right', 
        row_title_rot = 90, 
        row_gap = unit(0.25, 'cm'), 
        row_names_gp = gpar(fontsize=8,
                            col = c(rep("steelblue", 32), rep("indianred", 8)),
                            fontface = 'bold'),
        cluster_rows = T, 
        row_km = 4,
        row_dend_width = unit(3, 'cm'),
        cluster_columns = T, 
        show_column_dend = F, 
        column_split = factor(c(rep('PlasmidFinder', 28), 
                                rep('AST', 6),
                                rep("CTX-M-1", 4),
                                rep("CTX-M-9", 2),
                                rep("ß-lactam", 6),
                                rep("Agl", 8),
                                rep("Phe", 2),
                                rep("Tmt", 3),
                                rep("MLS", 3),
                                rep("Fcn", 1),
                                rep("Col", 2),
                                rep("Flq", 1),
                                rep("Sul", 3),
                                rep("Tet", 3),
                                rep('MLST', 9))), 
        column_title_gp = gpar(fontsize=8, fontface='bold'),
        column_gap = unit(0.5, 'cm'),
        heatmap_legend_param = list(title='', at = c(0,1),
                                    labels = c('Absent', 'Present')),
        rect_gp = gpar(col='white', lwd = 2))
dev.off()

save.image()

###########################################
# 210814 ordination plot based on resistome
###########################################
library(vegan)
library(ape)

# merge cordf.all5 with virulence factor
df.cor1 = cordf.all5 %>% rownames_to_column()

df.vf1 = df.count.vf %>% mutate(rowname = paste0(setting, '-', Sample)) %>% select(!c(setting, Sample))
df.vf1$rowname = gsub("JD-", "", df.vf1$rowname)

df.cor2 = inner_join(df.cor1, df.vf1, by='rowname') %>% column_to_rownames()

cor.dist = vegdist(cordf.all5,  method = "jaccard") # generate distance matrix
cor.pcoa = ape::pcoa(cor.dist) # generate principal component based on eigenvalues
barplot(cor.pcoa$values$Relative_eig[1:5]) # see variance explained by each component

# PCOA for resistome
cor.pcoa$vectors %>% data.frame %>% dplyr::select(Axis.1, Axis.2) %>% rownames_to_column() %>% 
  mutate(group = ifelse(grepl('community', rowname), 'Community', 'Clinical')) %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=group)) + geom_point() + theme_bw() +
  ggsci::scale_color_d3(name='Setting') + labs(x='Axis 1 (21.63% variance)', y='Axis 2 (15.23% variance)') 
  
ggsave('Fig5.tiff', units = 'in', height = 7, width=7, dpi=300, compression='lzw')

save.image()
# 210725 CTXM distribution of all WGS
name.clin

fig.ctxm = df.arg.tree2 %>% rownames_to_column() %>% mutate(setting = ifelse(rowname %in% name.clin, 'Clinical', 'Community')) %>% pivot_longer(contains('CTX')) %>% mutate(value = ifelse(value == 'Absent', 0, 1)) %>% group_by(name, setting) %>% summarise(n = sum(value)) %>% filter(n !=0) %>% ggplot(aes(x=reorder(name, n), y=n, fill=setting)) + geom_col() + facet_wrap(~setting, scales='free_x') + scale_y_continuous(breaks = seq(from = 0, to = 10, by = 1)) + theme_bw() + labs(x = 'CTX-M Variant', y = 'Count') + ggsci::scale_fill_d3() + theme(legend.position = 'none', axis.text.x = element_text(angle=45, hjust=1))

# 210725 MLST distribution of all WGS
fig.mlst = srst %>% select(ST) %>% filter(ST != 'NF*') %>% rownames_to_column() %>% mutate(setting = ifelse(rowname %in% name.clin, 'Clinical', 'Community')) %>% select(!rowname) %>% count(ST, setting) %>% ggplot(aes(x=reorder(ST,n), y=n, fill=ST)) + geom_col() + facet_wrap(~setting, scales='free_y')+ coord_flip() + theme_bw() + ggsci::scale_fill_igv() + theme(legend.position = 'none') + labs(x='ST', y='Count')

library(ggpubr)
ggarrange(fig.ctxm, fig.mlst, nrow=2, labels = 'auto')
ggsave('fig.ctxm.mlst.pdf', units='in', height=10,width=010)
ggsave('Fig2.tiff', units = 'in', height = 7, width=7, dpi=300, compression='lzw')
###############################
# risk factor of ESBL with meds
###############################
meds = read.csv('medication/qs_meds2.csv', header=T)
meds2 = meds %>% mutate(across(.cols = med_suppl_name, tolower),
                amlodipine = ifelse(med_suppl_name == 'amlodipine', 1, 0),
                simvastatin = ifelse(med_suppl_name == 'simvastatin', 1, 0),
                metformin = ifelse(med_suppl_name == 'metformin', 1, 0),
                meds_any = ifelse(med_suppl_name != '', 1, 0)) 

meds3 = df %>% mutate(id = paste0('JD-', stoolcode)) %>% select(id, hh, esbl)

meds4 = inner_join(meds2, meds3, by='id')

library(lme4)
library(boot)

mod1 = glmer(esbl ~ meds_any + (1|hh), data = meds4, family=binomial(link='logit'))
mod0 = glmer(esbl ~ 1 + (1|hh), data = meds4, family=binomial(link='logit'))

anova(mod1, mod0)
summary(mod1)

# comorbidities risk factor

cob = read.csv('medication/qs_meds.csv', header=T) %>% select(id, diabetes, cholesterol, bp, surgery_type) %>% 
  mutate(surgery = ifelse(surgery_type != '', 1, 0))
cob2 = inner_join(cob, meds3, by='id')

mod1 = glmer(esbl ~ surgery + (1|hh), data = cob2, family=binomial(link='logit'))
mod0 = glmer(esbl ~ 1 + (1|hh), data = cob2, family=binomial(link='logit'))

anova(mod1, mod0)

# CTXM association with CAZ

logit.caz = cordf.all5 %>% select(CAZ30, contains('CTX.M'))
logit.caz 
library(lme4)
library(boot)
library(performance)
library(see)
library(qqplotr)

model.caz = lm(CAZ30 ~ log(CTX.M.65 + 0.001), data = logit.caz) 
fig.modelcaz = check_model(model.caz) 
pdf(file = 'fig_model_caz.pdf', width=9, height=9)
fig.modelcaz
dev.off()
ggsave('fig_model_caz.pdf', units = 'in', height=9, width=9)
model.caz %>% summary

# ST131 typer

st131 = read.table('meta_analysis/st131typer.txt', sep='\t', header=T)
st131 %>% filter(Sample == 'JD-34.fasta')

####################################
# 211015 MANUSCRIPT REVISION PLOSONE
####################################
library(tidyverse)
library(corrplot)

df.count.arg # genotype list
df.ast.wgs # community ast phenotype

# isolate community isolates only
testdf = df.count.arg %>% filter(setting =='community') %>% select(-setting) 
phenotype_list = df.ast.wgs %>% rename(Sample = newcode) %>% inner_join(testdf, by='Sample') %>% select(1:13)
genotype_list = df.ast.wgs %>% rename(Sample = newcode) %>% inner_join(testdf, by='Sample') %>% select(-c(1:14))

library(corrplot)
testcor = df.ast.wgs %>% rename(Sample = newcode) %>% inner_join(testdf, by='Sample') %>% select(-Sample) %>% cor 
testcor.p = df.ast.wgs %>% rename(Sample = newcode) %>% inner_join(testdf, by='Sample') %>% select(-Sample) %>% cor.mtest()
testcor2 = testcor[-c(1:13),1:13]
testcor.p2 = testcor.p$p[-c(1:13),1:13]

# rough corrplot
corrplot(testcor2, method = 'circle')

# remove KZ30, IPM10, TZP110, CTX30 due to collinearity
# remove AMPH_Ecoli_Bla, AmpC2_Ecoli_Bla, ErmC_MLS, MrdA_Bla due to collinearity
testcor3 = testcor2[-c(1,7,17,24),-c(3,6,10,12)]
testcor.p3 = testcor.p2[-c(1,7,17,24),-c(3,6,10,12)]

rownames(testcor3) = rownames(testcor3) %>% str_split_fixed(pattern = '_', n=2) %>% data.frame %>% mutate(X2 = ifelse(.$X2 == 'AGly', 'Aminoglycoside', ifelse(
  .$X2 == 'Ecoli_Bla', 'Beta-lactam', ifelse(
    .$X2 == 'Bla', 'Beta-lactam', ifelse(
      .$X2 == 'Phe', 'Phenicols', ifelse(
        .$X2 == 'Tmt', 'Trimethoprim', ifelse(
          .$X2 == 'Fcyn', 'Fosfomycin', ifelse(
            .$X2 == 'MLS', 'Macrolides, Lincosamides, Streptogramines', ifelse(
              .$X2 == 'Flq', 'Fluoroquinolone', ifelse(
                .$X2 == 'Sul', 'Sulfonamide', ifelse(
                  .$X2 == 'Tet', 'Tetracycline', 'Colistin'))))))))))) %>% mutate(
                    revised = paste0(X2, '_', X1)) %>% pull(revised)

rownames(testcor.p3) = rownames(testcor3)

sorted.name = rownames(testcor3) %>% sort
testcor3[sorted.name,]

pdf(file = 'draft4_submission/revision1/fig_ast_arg_corr.pdf', height=6, width=8)
corrplot(testcor3[sorted.name,], 
         p.mat = testcor.p3[sorted.name,], 
         method = 'color',
         insig = 'label_sig',
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = 'black',
         pch.cex = 1)
dev.off()

# no of significant association per phenotype
testcor.p3 %>% data.frame %>% rownames_to_column() %>% pivot_longer(!rowname) %>% filter(value < 0.05) %>% count(name) %>% arrange(-n)
# no of significant association per phenotype
testcor.p3 %>% data.frame %>% rownames_to_column() %>% pivot_longer(!rowname) %>% filter(value < 0.05) %>% count(rowname) 

testcor.p3 %>% data.frame %>% rownames_to_column() %>% pivot_longer(!rowname) %>% filter(value < 0.05) %>% filter(name == 'CIP5')

# plos one revision pointfinder

point = read.table('df_pointfinder.tsv', header=T, sep = '\t')
name.parc = point %>% filter(grepl('parC', .$Mutation)) %>% pull(sample)
name.pare = point %>% filter(grepl('parE', .$Mutation)) %>% pull(sample)
name.gyra = point %>% filter(grepl('gyr', .$Mutation)) %>% pull(sample)

name.qnrs = srst.all %>% filter(grepl('Qnr', .$name)) %>% filter(value2 == 'present') %>% pull(Sample)

# annotate ast with qrdr and pmqr

fluor.clin = df.ast.all %>% filter(Sample %in% name.clin) %>% select(Sample,CIP5) %>% mutate(gyrA = ifelse(.$Sample %in% name.gyra, 1, 0),
                                                                                parC = ifelse(.$Sample %in% name.parc, 1, 0),
                                                                                parE = ifelse(.$Sample %in% name.pare, 1, 0),
                                                                                qnrS = ifelse(.$Sample %in% name.qnrs, 1, 0))

fluor.com = df.ast.wgs2 %>% rownames_to_column('Sample') %>% select(Sample, CIP5, NA30) %>% mutate(gyrA = ifelse(.$Sample %in% name.gyra, 1, 0),
                                                                                       parC = ifelse(.$Sample %in% name.parc, 1, 0),
                                                                                       parE = ifelse(.$Sample %in% name.pare, 1, 0),
                                                                                       qnrS = ifelse(.$Sample %in% name.qnrs, 1, 0))


# annotate with mlst
name.131 = c(srst.st.com %>% filter(ST == 131) %>% pull(newcode),
             srst %>% rownames_to_column('sample') %>% filter(sample %in% name.clin) %>% filter(ST == 131) %>% pull(sample)
)

fluor.clin2 = fluor.clin %>% pivot_longer(!Sample) %>% mutate(lab = ifelse(name == 'CIP5' | name == 'NA30', 'AST', ifelse(
  name == 'qnrS', 'PMQR', 'QRDR'))) 

fluor.com2 = fluor.com %>% pivot_longer(!Sample) %>% mutate(lab = ifelse(name == 'CIP5' | name == 'NA30', 'AST', ifelse(
  name == 'qnrS', 'PMQR', 'QRDR'))) 

library(ggtext)
library(glue)
library(ggh4x)

bind_rows(fluor.clin2, fluor.com2) %>% 
  mutate(Setting = ifelse(Sample %in% name.clin, 'Clinical', 'Community')) %>%
  mutate(name = ifelse(lab == 'AST', name, paste0('*', name, '*'))) %>% 
  ggplot(aes(x = Sample, y = name, fill = as.factor(value))) + 
                                                  geom_tile() + 
                                                  scale_fill_manual(name = 'Positivity', 
                                                                    labels = c('Negative', 'Positive'), 
                                                                    values = c('steelblue', 'orange')) + coord_flip() + 
  facet_grid2(rows = vars(Setting), cols = vars(lab), scales='free') +
  theme_bw() +
  labs(y = '', x='Samples') +
  theme(axis.text.x = element_markdown(),
        axis.text.y = element_markdown(size = 6))
ggsave('draft4_submission/revision1/fig_flq.pdf', units = 'in', height=6,width=8)


# no of isolates with specific PMQR or QRDR mutation
bind_rows(fluor.clin2, fluor.com2) %>% 
  mutate(Setting = ifelse(Sample %in% name.clin, 'Clinical', 'Community')) %>%
  filter(lab == 'QRDR',
         value == 1,
         name =='parE')

# phenotypic-genotypic correlation?
df.flq = bind_rows(
fluor.com %>% mutate(pheno = CIP5 + NA30,
                     geno = gyrA + parC + parE + qnrS) %>% 
  select(Sample, pheno, gyrA, parC, parE, qnrS)
,
fluor.clin %>% mutate(pheno = CIP5,
                 geno = gyrA + parC + parE + qnrS) %>%
  select(Sample, pheno, gyrA, parC, parE, qnrS)
) %>% mutate(pheno = ifelse(pheno >0, 1, 0))

library(corrplot)
flq.cor = df.flq %>% column_to_rownames('Sample') %>% cor
flq.sig = cor.mtest(flq.cor)
corrplot(flq.cor, p.mat = flq.sig$p)

# no isolates with mutation
qrdr.name = point %>% pull(sample) %>% unique

# compare with ast output for fluoroquinolone
cip.name = df.ast.all %>% filter(CIP5 == 1) %>% pull(Sample)
cip.name %in% qrdr.name %>% table

#exPEC marker based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC161843/
expec.name = srst.vf3 %>% 
  select(!setting) %>% 
  rownames_to_column() %>% select(c(rowname,
                                    contains('papA'), 
                                    contains('papC'),
                                    contains('iut'), 
                                    contains('kps'), 
                                    contains('afa'), 
                                    contains('dra'), 
                                    contains('sfa'), 
                                    contains('foc'))) %>%
  mutate(expec.pap = ifelse(papA + papC > 0, 1, 0),
         expec.kps = ifelse(kpsC + kpsD + kpsE + kpsF + kpsM + kpsS + kpsU > 0, 1, 0),
         expec.iut = ifelse(iutA > 0, 1, 0),
         expec.afa.dra = ifelse(afaF.III + draA + draB + draC + draD + draP > 0, 1, 0)) %>%
  select(c(rowname, contains('expec'))) %>%
  column_to_rownames('rowname') %>%
  mutate(sum = rowSums(.),
         expec.st = ifelse(sum >= 2, 1, 0)) %>%
  rownames_to_column() %>%
  filter(expec.st == 1) %>% pull(rowname)

# UPEC marker based on https://journals.asm.org/doi/full/10.1128/IAI.00752-12 
# and https://www.frontiersin.org/articles/10.3389/fmicb.2017.00146/full#B35

upec.name = srst.vf3 %>%
  select(!setting) %>% 
  rownames_to_column() %>% select(c(rowname,
                                    contains('fyuA'), 
                                    contains('yfc'),
                                    contains('chuA'), 
                                    contains('vat'), 
                                    contains('foc'), 
                                    contains('pap'), 
                                    contains('sfa'), 
                                    contains('cnf'))) %>%
  mutate(upec.fyua = fyuA,
         upec.chuA = chuA,
         upec.vat = vat,
         upec.pap = ifelse(papA + papB + papC + papD + papE + papF + papG + papH + papI + papJ + papK + papX > 0, 1, 0),
         upec.cnf = cnf1) %>% 
  select(rowname, contains('upec')) %>%
  mutate(upec.st = upec.fyua + upec.chuA + upec.vat + upec.pap + upec.cnf) %>%
  select(rowname, upec.st) %>%
  filter(upec.st > 0) %>%
  pull(rowname)
  
# EAEC marker based on https://www.frontiersin.org/articles/10.3389/fmicb.2017.00146/full#B35

eaec.name = srst.vf3 %>%
  select(!setting) %>% 
  rownames_to_column() %>% select(c(rowname,
                                    contains('aatA'), 
                                    contains('agg'),
                                    contains('aap'))) %>%
  mutate(eaec.st = aatA + aap.aspU) %>%
  filter(eaec.st == 1) %>%
  pull(rowname)
                                    
# EHEC marker based on https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-11-142#Sec7
# no STEC = stx gene = shiga-toxin
# no atypical EPEC = eae gene - eaeX = invasin/intimin homologue

epec.atypical.name = srst.vf3 %>%
  select(!setting) %>% 
  rownames_to_column() %>% select(c(rowname,
                                    contains('eae'))) %>%
  rename(epec.atypical.st = air.eaeX) %>% 
  filter(epec.atypical.st == 1) %>%
  pull(rowname)
  
# compile pathogroups
library(ggtext)
srst.vf3 %>% rownames_to_column() %>%
  select(rowname, setting) %>%
  mutate(expec = ifelse(rowname %in% expec.name, 1, 0),
         upec = ifelse(rowname %in% upec.name, 1, 0),
         eaec = ifelse(rowname %in% eaec.name, 1, 0),
         epec.atypical = ifelse(rowname %in% epec.atypical.name,1 ,0)) %>% 
  select(!setting) %>%
  pivot_longer(!rowname) %>% 
  mutate(setting = ifelse(rowname %in% name.clin, 'Clinical', 'Community')) %>%
  ggplot(aes(y=rowname, x=name, fill=as.factor(value))) + geom_tile() +
  scale_fill_manual(name = '',
                    labels = c('Negative', 'Positive'),
                    values = c('steelblue', 'orange')) +
  scale_x_discrete(limits = c('expec', 'upec', 'eaec', 'epec.atypical'),
                   labels = c('ExPEC', 'UPEC', 'EAEC', 'EPEC-<br>atypical')) +
  facet_wrap(~setting, scales = 'free') +
  theme_bw() +
  theme(axis.text.x = element_markdown()) +
  labs(y = '', x='')

ggsave('draft4_submission/revision1/fig_pathogroup.pdf', units='in', height=6,width=8)
  
  
srst.vf3 %>% rownames_to_column() %>%
  select(rowname, setting) %>%
  mutate(expec = ifelse(rowname %in% expec.name, 1, 0),
         upec = ifelse(rowname %in% upec.name, 1, 0),
         eaec = ifelse(rowname %in% eaec.name, 1, 0),
         epec.atypical = ifelse(rowname %in% epec.atypical.name,1 ,0)) %>% filter(upec == 1) %>% pull(setting) %>% table
