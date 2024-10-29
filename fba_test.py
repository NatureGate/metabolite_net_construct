from cobra.io import  load_model
import cobra
model = cobra.io.read_sbml_model('3sl/iECD_1391.xml')
reaction = model.reactions.get_by_id('23PDE2pp')
reaction.notes['gibbs'] = 33.3
print(model.reactions.get_by_id('23PDE2pp').annotation)




