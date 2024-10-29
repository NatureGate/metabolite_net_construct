# Stoichiometry
# BiGG ID
# -1.0 acgam1p_c/N-Acetyl-D-glucosamine 1-phosphate
# -1.0 h_c
# 1.0 ppi_c
# 1.0 uacgam_c  /UDP-N-acetyl-D-glucosamine
# -1.0 utp_c
# Name
# N-Acetyl-D-glucosamine 1-phosphate
# H+
# Diphosphate
# UDP-N-acetyl-D-glucosamine
# UTP C9H11N2015P3

import cobra
from cobra import Metabolite, Reaction
from cobra.io import write_sbml_model
from cobra.manipulation.delete import delete_model_genes

delete_genes = ['nagB', 'nagA', 'nanE', 'nanK', 'nanT', 'nanA', 'lacZ']
delete_gene_numbers = ['b0678', 'b0677', 'b3223', 'b3222', 'b3224', 'b3225', 'b0344']
bl21_model = cobra.io.read_sbml_model('iML1515.xml')
delete_model_genes(bl21_model, [delete_gene_numbers])
cobra.manipulation.remove_genes(bl21_model, ['b0678', 'b0677', 'b3223', 'b3222', 'b3224', 'b3225', 'b0344'])
UDP_Glc_NAc = bl21_model.metabolites.get_by_id('uacgam_c')
# if ''UDP_Glc_NAc.annotation.keys()
lcts_c = bl21_model.metabolites.get_by_id('lcts_c')
udpg_c = bl21_model.metabolites.get_by_id('udpg_c')
ex_lcts_e_reaction = bl21_model.reactions.get_by_id('EX_lcts_e')
ex_lcts_e_reaction.lower_bound = -10
########################################################################
gam6p_c = bl21_model.metabolites.get_by_id('gam6p_c')
h2o_c = bl21_model.metabolites.get_by_id('h2o_c')
f6p_c = bl21_model.metabolites.get_by_id('f6p_c')
nh4_c = bl21_model.metabolites.get_by_id('nh4_c')
# G6PDA_reaction = Reaction('G6PDA_reverse')
# G6PDA_reaction.name = 'G6PDA_reverse'
# G6PDA_reaction.lower_bound = 0
# G6PDA_reaction.upper_bound = 1000
# G6PDA_reaction.add_metabolites({
#     gam6p_c: 1.0,
#     h2o_c: 1.0,
#     f6p_c: -1.0,
#     nh4_c: -1.0
# })
# bl21_model.add_reactions([G6PDA_reaction])
########################################################################
# h2o_c + uacgam_c → h_c + udp_c + acmana_c
udp_glc_nac = bl21_model.metabolites.get_by_id('uacgam_c')  # UDP-Glc-NAc
man_nac = bl21_model.metabolites.get_by_id('acmana_c')  # man-NAc
neu5ac = bl21_model.metabolites.get_by_id('acnam_c')  # Neu5Ac
udp = bl21_model.metabolites.get_by_id('udp_c')
h = bl21_model.metabolites.get_by_id('h_c')

# h2o_c + uacgam_c → h_c + udp_c + acmana_c
mannac_reaction = cobra.Reaction("R_ManNAc")
mannac_reaction.name = "mannac_syn"
mannac_reaction.lower_bound = 0.0  #
mannac_reaction.upper_bound = 1000.0

# 定义反应式
metabolites = {
    h2o_c: -1.0,  # 底物
    udp_glc_nac: -1.0,  #
    man_nac: 1.0,
    udp: 1.0,
    h: 1.0
}

mannac_reaction.add_metabolites(metabolites)
# 添加新的基因
neuc_gene = cobra.Gene("neuC")
neuc_gene.name = "neuC"
bl21_model.genes.append(neuc_gene)

# 设置基因规则
mannac_reaction.gene_reaction_rule = "neuC"
bl21_model.add_reactions([mannac_reaction])
############################################################################################
# phosphoenolpyruvate + N-acetyl-D-mannosamine + H2O → N-acetyl-β-neuraminate + phosphate
#        PEP                        MANAc                              Neu5Ac
pep = bl21_model.metabolites.get_by_id('pep_c')  # PEP
pi = bl21_model.metabolites.get_by_id('pi_c')  # HO4P
h2o = bl21_model.metabolites.get_by_id('h2o_c')  # HO4P
neu5ac_reaction = cobra.Reaction("R_neu5ac")

neu5ac_reaction.name = "neu5ac_syn"
neu5ac_reaction.lower_bound = 0.0  #
neu5ac_reaction.upper_bound = 1000.0

# 定义反应式
neu5ac_metabolites = {
    man_nac: -1.0,  # 底物
    pep: -1.0,  #
    h2o: -1.0,
    neu5ac: 1.0,
    pi: 1.0
}
neu5ac_reaction.add_metabolites(neu5ac_metabolites)
# 添加新的基因
neuc_gene = cobra.Gene("neuB")
neuc_gene.name = "neuB"
bl21_model.genes.append(neuc_gene)

# 设置基因规则
neu5ac_reaction.gene_reaction_rule = "neuB"
bl21_model.add_reactions([neu5ac_reaction])

##########################################################################################

# CTP + N - Acetylneuraminate <= > Diphosphate + CMP - N - acetylneuraminate
# https://www.genome.jp/entry/R01117
ctp = bl21_model.metabolites.get_by_id('ctp_c')
ppi = bl21_model.metabolites.get_by_id('ppi_c')
# h2o = bl21_model.metabolites.get_by_id('cmpacneun_c')
cmp_neu5ac = Metabolite(
    'cmpacneun',
    name='cmpneu5ac',
    formula='C20H31N4O16P',
    compartment='c',
    charge=0
)
# 定义反应式
cmpneu5ac_metabolites = {
    ctp: -1.0,  # 底物
    neu5ac: -1.0,  #
    ppi: 1.0,
    cmp_neu5ac: 1.0
}
cmpneu5ac_reaction = cobra.Reaction("R_cmpneu5ac")
cmpneu5ac_reaction.name = "cmpneu5ac_syn"
cmpneu5ac_reaction.lower_bound = 0.0  #
cmpneu5ac_reaction.upper_bound = 1000.0
cmpneu5ac_reaction.add_metabolites(cmpneu5ac_metabolites)
bl21_model.add_reactions([cmpneu5ac_reaction])
# 添加新的基因
cmpneuc_gene = cobra.Gene("neuA")
cmpneuc_gene.name = "neuA"
bl21_model.genes.append(cmpneuc_gene)
# 设置基因规则
neu5ac_reaction.gene_reaction_rule = "neuA"
##############################################################################################
lcts_c = bl21_model.metabolites.get_by_id('lcts_c')
sl3 = Metabolite(
    'sl3',
    name='sl3',
    formula='C23H39NO19',
    compartment='c',
    charge=0
)
sl3_metabolites = {
    cmp_neu5ac: -1.0,  # 底物
    lcts_c: -1.0,  #
    sl3: 1.0,
    h2o: 1.0
}

sl3_reaction = cobra.Reaction("R_sl3")
sl3_reaction.name = "sl3_syn"
sl3_reaction.lower_bound = 0.0  #
sl3_reaction.upper_bound = 1000.0
sl3_reaction.add_metabolites(sl3_metabolites)
bl21_model.add_reactions([sl3_reaction])

# 添加新的基因
sl3_gene = cobra.Gene("vs16")
sl3_gene.name = "vs16"
bl21_model.genes.append(sl3_gene)
# 设置基因规则
sl3_reaction.gene_reaction_rule = "vs16"
# 添加demand反应
# bl21_model.add_boundary(sl3, type="exchange")
bl21_model.add_boundary(sl3, type="demand")
bl21_model.objective = "R_sl3"
print(bl21_model.objective)
print(bl21_model.optimize())

# to_del = bl21_model.genes.get_by_id('b3223')
# print(to_del)

write_sbml_model(bl21_model, "sl3_withnot_knockout.xml")
