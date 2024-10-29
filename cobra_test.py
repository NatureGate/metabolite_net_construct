import networkx as nx
import cobra 
import pandas as pd
import time
import sys
use_model_file="data/external/iML1515.xml" 
model = cobra.io.read_sbml_model(use_model_file)
cobra_config = cobra.Configuration()
##reaction bounds
solver = model.solver
# for reaction in model.reactions:
#     for (key,value) in reaction.metabolites.items():
#         print(key,value)
#         model.metabolites
#         break
#     break
from cobra import Model,Reaction,Metabolite
#创建一个model
model = Model('example_model')
#创建一个reaction
reaction = Reaction('R_3OAS140')
reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000
##创建代谢物
ACP_c = Metabolite(
    'ACP_c',
    formula='C11H21N2O7PRS',
    name='acyl-carrier-protein',
    compartment='c')
omrsACP_c = Metabolite(
    'M3omrsACP_c',
    formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment='c')
co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
malACP_c = Metabolite(
    'malACP_c',
    formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein',
    compartment='c')
h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
ddcaACP_c = Metabolite(
    'ddcaACP_c',
    formula='C23H43N2O8PRS',
    name='Dodecanoyl-ACP-n-C120ACP',
    compartment='c')
reaction.add_metabolites({
    malACP_c: -1.0,
    h_c: -1.0,
    ddcaACP_c: -1.0,
    co2_c: 1.0,
    ACP_c: 1.0,
    omrsACP_c: 1.0
})
model.add_reactions([reaction])
reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
reaction.genes

# The objects have been added to the model
print(f'{len(model.reactions)} reactions')
print(f'{len(model.metabolites)} metabolites')
print(f'{len(model.genes)} genes')
model.add_metabolites([Metabolite('glycogen_c',name='glycogen',compartment='c'),Metabolite('co2_e',name='CO2',compartment='e'),])
model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("glycogen_c"), type="sink")
print("exchanges", model.exchanges)
print("sinks", model.sinks)
print("demands", model.demands)

# boundary reactions
print(model.boundary)
set(model.reactions) - set(model.boundary)




from cobra.io import load_model
model = load_model("textbook")


