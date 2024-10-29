import cobra
from cobra import Metabolite, Reaction
from cobra.io import write_sbml_model
from cobra.manipulation.delete import delete_model_genes
import cobra
from cobra.io import read_sbml_model

# delete_genes = ['nagB','nagA','nanE','nanK','nanT','nanA','lacZ']
# delete_gene_numbers = ['b0678','b0677','b3223','b3222','b3224','b3225','b0344']
# bl21_model = cobra.io.read_sbml_model('iML1515.xml')
# # delete_model_genes(bl21_model, [delete_gene_numbers])
# # cobra.manipulation.knock_out_model_genes(bl21_model, ['b0678', 'b0677', 'b3223', 'b3222', 'b3224', 'b3225', 'b0344'])
# UDP_Glc_NAc = bl21_model.metabolites.get_by_id('uacgam_c')
# b3222_gene = bl21_model.genes.get_by_id('b3222')
# b3222_reactions = b3222_gene.reactions
#
# b3224_gene = bl21_model.genes.get_by_id('b3224')
# b3224_reactions = b3224_gene.reactions
# print(b3224_reactions)

model = read_sbml_model('sl3.xml')
# org_sl3_obj_val = model.optimize().objective_value
# model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
# org_biomass_obj_val = model.optimize().objective_value
# genes = model.genes
#
# for gene in genes:
#     gene_id = gene.id
#     new_model = model
#     cobra.manipulation.remove_genes(new_model, [gene_id])
#     biomass_obj_val = new_model.optimize().objective_value
#     new_model.objective = 'R_sl3'
#     obj_val = new_model.optimize().objective_value
#     if obj_val > org_sl3_obj_val:
#         print(gene_id)
#         if biomass_obj_val / org_biomass_obj_val >= 0.8:
#
#             print(gene_id)
# reaction_del = 'CYTBO3_4pp'
# Genes:
# b0429 (cyoD)
# b0432 (cyoA)
# b0431 (cyoB)
# b0430 (cyoC)
# gene = ['b0429', 'b0430', 'b0431', 'b0432']
# model.objective = "R_sl3"
# print(model.objective)
# print(model.optimize())
# cobra.manipulation.remove_genes(model, ['b1278'])
# print(model.objective)
# print(model.optimize())
model = read_sbml_model('sl3.xml')
model.objective = 'DM_sl3'
obj_val = model.optimize().objective_value
print(obj_val)