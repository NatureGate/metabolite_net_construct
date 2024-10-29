import cobra
model = cobra.io.read_sbml_model('lnfp.xml')
# model = cobra.io.mat.load_matlab_model('lnfp.mat')
# reaction = model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M')
# print(model.solver)
for reaction in model.reactions:
    if reaction.id.__contains__('MASS'):
        print(reaction.id)

for metabolite in model.metabolites:
    mt = metabolite