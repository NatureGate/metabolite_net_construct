import cobra
from cobra import Metabolite, Reaction
from cobra.io import write_sbml_model, save_matlab_model
import cobra.manipulation
from cobra.io import load_json_model, save_json_model, load_matlab_model

bl21_model = cobra.io.read_sbml_model('iECD_1391.xml')
gdp_c = bl21_model.metabolites.get_by_id('gdp_c')
gdp_l_fuc = bl21_model.metabolites.get_by_id('gdpfuc_c')
lcts_c = bl21_model.metabolites.get_by_id('lcts_c')

ex_lcts_e_reaction = bl21_model.reactions.get_by_id('EX_lcts_e')
ex_lcts_e_reaction.lower_bound =-10
# 转运乳糖，2FL还有一个底物是GDP-L-Fuc，GDP-L-Fuc跟Lactose合成2FL，然后生成一个GDP
reaction = Reaction('fl2')
reaction.name = '2′-Fucosyllactose synthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000

# 添加2FL细胞内
FL2_c = Metabolite(
    'FL2_c',
    formula='C18H32O15',
    name='2-fucosyllactose c',
    compartment='c')
# 添加2FL细胞外
FL2_e = Metabolite(
    'FL2_e',
    formula='C18H32O15',
    name='2-fucosyllactose e',
    compartment='e')

reaction.add_metabolites({
    lcts_c: -1.0,
    gdp_l_fuc: -1.0,
    FL2_c: 1.0,
    gdp_c: 1.0
})

metabolite_map = reaction.metabolites
bl21_model.add_reactions([reaction])
# 添加2FL
bl21_model.add_metabolites([FL2_c,FL2_e])
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("FL2_e"), type="exchange")
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("FL2_c"), type="demand")
reaction.gene_reaction_rule = 'alpha12fuct'





##构建LNT

UDP_Glc_NAc = bl21_model.metabolites.get_by_id('uacgam_c')
lcts_c = bl21_model.metabolites.get_by_id('lcts_c')
udpg_c = bl21_model.metabolites.get_by_id('udpg_c')
ex_lcts_e_reaction = bl21_model.reactions.get_by_id('EX_lcts_e')
ex_lcts_e_reaction.lower_bound = -10

##3糖合成
LNTriII_reaction = Reaction('LNTriII')
LNTriII_reaction.name = 'Lacto‑N‑Tetraose'
LNTriII_reaction.lower_bound = 0
LNTriII_reaction.upper_bound = 1000

# 添加LNTriII细胞外
LNTriII_e = Metabolite(
    'LNTriII_e',
    formula='C20H35NO16',
    name='Lacto-N-triose II e',
    compartment='e')
# 添加LNTriII到细胞内
LNTriII_c = Metabolite(
    'LNTriII_c',
    formula='C20H35NO16',
    name='Lacto-N-triose II c',
    compartment='c')
# LNTriII 合成反应
LNTriII_reaction.add_metabolites({
    lcts_c: -1.0,
    UDP_Glc_NAc: -1.0,
    LNTriII_c: 1.0,

})
bl21_model.add_reactions([LNTriII_reaction])

# 添加3糖的边界反应
bl21_model.add_metabolites([LNTriII_e, LNTriII_c])
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("LNTriII_e"), type="exchange")
bl21_model.add_boundary(LNTriII_c, type="demand")

LNTriII_reaction.gene_reaction_rule = 'lgtA'

#####################################################################################################
#
# 先添加UDP-Gal反应 UDP-Glc<-->UDP-Gal
######################################################################################################
udpg_c = bl21_model.metabolites.get_by_id('udpg_c')
udpgal_c = Metabolite(
    'udpgal',
    name='UDPgalactose',
    formula='C15H22N2O17P2',
    compartment='c',
    charge=-2
)

UDPG4E_reaction = Reaction('UDPG4E')
UDPG4E_reaction.name = 'UDPglucose 4-epimerase'
UDPG4E_reaction.lower_bound = -1000
UDPG4E_reaction.upper_bound = 1000
UDPG4E_reaction.add_metabolites({
    udpg_c: -1.0,
    udpgal_c: 1.0,
})
UDPG4E_reaction.gene_reaction_rule = 'galE'
bl21_model.add_reactions([UDPG4E_reaction])
#####################################################################################################
#
# 添加LNT合成反应
######################################################################################################
lnt_c = Metabolite(
    'lnt_c',
    name='lnt c',
    formula='C26H45NO21',
    compartment='c',

)

lnt_e = Metabolite(
    'lnt_e',
    name='lnt e',
    formula='C26H45NO21',
    compartment='e',
)
##LNT合成
lnt_reaction = Reaction('LNT')
lnt_reaction.name = 'Lacto‑N‑Tetraose'
lnt_reaction.lower_bound = 0
lnt_reaction.upper_bound = 1000
lnt_reaction.add_metabolites({
    udpgal_c: -1,
    LNTriII_c: -1,
    lnt_c: 1
})
lnt_reaction.gene_reaction_rule = 'beta13galT'
bl21_model.add_metabolites([lnt_e, lnt_c])
bl21_model.add_reactions([lnt_reaction])
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("lnt_c"), type="demand")
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("lnt_e"), type="exchange")
# bl21_model.objective = "LNTriII"
# bl21_model.objective = "fl2"
# print(bl21_model.genes.get_by_id('ECD_00298'))
# cobra.manipulation.knock_out_model_genes(model=bl21_model,gene_list=['b2047','b0344','b3916'])
solver = bl21_model.optimize()
# write_sbml_model(bl21_model, "LNT.xml")
print(solver)


# LNFP生成

# 添加LNTriII细胞外
lnfp1_e = Metabolite(
    'LNFP1',
    formula='C32H55NO25',
    name='LNFP1',
    compartment='e')
# 添加lnfp到细胞内
lnfp1_c = Metabolite(
    'LNFP1_c',
    formula='C32H55NO25',
    name='LNFP1',
    compartment='c')

bl21_model.add_boundary(lnfp1_e, type="exchange")
bl21_model.add_boundary(lnfp1_c, type="demand")

lnfp_reaction = Reaction('lnfp')
lnfp_reaction.name = 'LNFP'
lnfp_reaction.lower_bound = 0
lnfp_reaction.upper_bound = 1000
lnfp_reaction.add_metabolites({
    gdp_l_fuc: -1,
    lnt_c: -1,
    lnfp1_c:1,
    gdp_c:1
})
bl21_model.add_reactions([lnfp_reaction])
lnfp_reaction.gene_reaction_rule = 'alpha12fuct'
bl21_model.objective = "lnfp"
print(bl21_model.objective)
print(bl21_model.optimize())
# write_sbml_model(bl21_model, "lnfp.xml")
#
save_matlab_model(bl21_model, "lnfp.mat")
