import cobra
from cobra import Metabolite, Reaction
from cobra.io import write_sbml_model, save_matlab_model,load_matlab_model
import cobra.manipulation

mini_mat_path = "lnfp.mat"
model = load_matlab_model(mini_mat_path)
print()
