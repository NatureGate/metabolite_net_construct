import cobra
from cobra.io import  save_matlab_model
model = cobra.io.read_sbml_model('sl3.xml')
save_matlab_model(model,'sl3.mat')