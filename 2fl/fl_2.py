import cobra
from cobra import Metabolite, Reaction
from cobra.io import write_sbml_model
##1515构建2fl
model = cobra.io.read_sbml_model("../data/external/iML1515.xml")
gdp_c = model.metabolites.get_by_id('gdp_c')
gdp_l_fuc = model.metabolites.get_by_id('gdpfuc_c')
lcts_c = model.metabolites.get_by_id('lcts_c')


ex_lcts_e_reaction = model.reactions.get_by_id('EX_lcts_e')
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
model.add_reactions([reaction])
# 添加2FL
model.add_metabolites([FL2_c,FL2_e])
model.add_boundary(model.metabolites.get_by_id("FL2_e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("FL2_c"), type="demand")
#model.add_boundary(model.metabolites.get_by_id("gdpfuc_c"), type="demand")
#model.add_boundary(model.metabolites.get_by_id("lcts_c"), type="demand")

# 目标修改为DM_gdpfuc_c时，solution是3.275，说明能到gdbfuc这一步
#model.objective = "DM_gdpfuc_c"
# 目标修改为DM_lcts_c
#model.objective = "DM_lcts_c"

# write_sbml_model(model, "FL2.xml")
#model.objective = "fl2"
#model.objective = "MAN6PI"
# model.objective = "LCTStpp"
solver = model.optimize()
#model.reactions.get_by_id("fl2").upper_bound = 1000.
print(solver)



