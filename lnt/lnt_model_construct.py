import cobra
from cobra import Metabolite, Reaction
from cobra.io import write_sbml_model

##构建LNT
bl21_model = cobra.io.read_sbml_model('iECD_1391.xml')
UDP_Glc_NAc = bl21_model.metabolites.get_by_id('uacgam_c')
lcts_c = bl21_model.metabolites.get_by_id('lcts_c')
udpg_c = bl21_model.metabolites.get_by_id('udpg_c')
ex_lcts_e_reaction = bl21_model.reactions.get_by_id('EX_lcts_e')
ex_lcts_e_reaction.lower_bound = -10

##3糖合成
reaction = Reaction('LNTriII')
reaction.name = 'Lacto‑N‑Tetraose'
reaction.lower_bound = 0
reaction.upper_bound = 1000

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
reaction.add_metabolites({
    lcts_c: -1.0,
    UDP_Glc_NAc: -1.0,
    LNTriII_c: 1.0,

})
bl21_model.add_reactions([reaction])

# 添加3糖的边界反应
bl21_model.add_metabolites([LNTriII_e, LNTriII_c])
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("LNTriII_e"), type="exchange")

reaction.gene_reaction_rule = 'lgtA'

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
UDPG4E_reaction.gene_reaction_rule = 'beta13galT'
bl21_model.add_metabolites([lnt_e, lnt_c])
bl21_model.add_reactions([lnt_reaction])
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("lnt_c"), type="demand")
bl21_model.add_boundary(bl21_model.metabolites.get_by_id("lnt_e"), type="exchange")
# write_sbml_model(bl21_model, "LNT.xml")
bl21_model.objective = "LNT"
solver = bl21_model.optimize()


print(solver)
