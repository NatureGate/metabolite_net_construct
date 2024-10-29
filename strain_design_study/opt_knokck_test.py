import straindesign as sd
import cobra

ecc = cobra.io.load_model('e_coli_core')
# Create copy of model to which pathway will be added
ecc_14bdo = ecc.copy()

# Add metabolites to model
ecc_14bdo.add_metabolites([cobra.Metabolite('sucsal_c'),  # Succinic semialdehyde
                           cobra.Metabolite('4hb_c'),  # 4-Hydroxybutanoate
                           cobra.Metabolite('4hbcoa_c'),  # 4-Hydroxybutyryl-CoA
                           cobra.Metabolite('4hbal_c'),  # 4-Hydroxybutanal
                           cobra.Metabolite('14bdo_c'),  # Butane-1,4-diol (cytopl.)
                           cobra.Metabolite('14bdo_p'),  # Butane-1,4-diol (peripl.)
                           cobra.Metabolite('14bdo_e')  # Butane-1,4-diol (extrac.)
                           ])

# Create reactions
SSCOARx = cobra.Reaction('SSCOARx')
AKGDC = cobra.Reaction('AKGDC')
HBD = cobra.Reaction('4HBD')
HBCT = cobra.Reaction('4HBCT')
HBDH = cobra.Reaction('4HBDH')
HBDx = cobra.Reaction('4HBDx')
BDOtpp = cobra.Reaction('14BDOtpp')
BDOtex = cobra.Reaction('14BDOtex')
EX_14bdo_e = cobra.Reaction('EX_14bdo_e')

# Add reactions to model
ecc_14bdo.add_reactions([SSCOARx,
                         AKGDC,
                         HBD,
                         HBCT,
                         HBDH,
                         HBDx,
                         BDOtpp,
                         BDOtex,
                         EX_14bdo_e])

# Define reaction equations
SSCOARx.reaction = '1 h_c + 1 nadph_c + 1 succoa_c -> 1 coa_c + 1 nadp_c + 1 sucsal_c'
AKGDC.reaction = '1 akg_c + 1 h_c -> 1 co2_c + 1 sucsal_c'
HBD.reaction = '1 h_c + 1 nadh_c + 1 sucsal_c  -> 1 4hb_c + 1 nad_c'
HBCT.reaction = '1 4hb_c + 1 accoa_c            -> 1 4hbcoa_c + 1 ac_c'
HBDH.reaction = '1 4hbcoa_c + 1 h_c + 1 nadh_c  -> 1 4hbal_c + 1 coa_c + 1 nad_c'
HBDx.reaction = '1 4hbal_c + 1 h_c + 1 nadh_c   -> 1 14bdo_c + 1 nad_c'
BDOtpp.reaction = '1 14bdo_c                      -> 1 14bdo_p'
BDOtex.reaction = '1 14bdo_p                      -> 1 14bdo_e'
EX_14bdo_e.reaction = '1 14bdo_e                      ->'

# Verify that pathway is operational
sol = sd.fba(ecc_14bdo, obj='EX_14bdo_e', obj_sense='max')
print(f"Maximum possible 1,4-BDO synthesis rate: {sol.objective_value}.")

module_optknock = sd.SDModule(ecc_14bdo, sd.names.OPTKNOCK,
                              inner_objective='BIOMASS_Ecoli_core_w_GAM',
                              outer_objective='EX_14bdo_e',
                              constraints='BIOMASS_Ecoli_core_w_GAM >= 0.5')

import logging

logging.basicConfig(level=logging.INFO)
## Compute strain designs
# allow all gene knockouts except for spontanuos
gko_cost = {g.name: 1 for g in ecc_14bdo.genes}
gko_cost.pop('s0001')
# possible knockout of O2
ko_cost = {'EX_o2_e': 1}
# addition candidates
ki_cost = {'AKGDC': 1, 'SSCOARx': 1}  # AKGDC was added in example 1.c)

sols = sd.compute_strain_designs(ecc_14bdo,
                                 sd_modules=module_optknock,
                                 max_solutions=1,
                                 max_cost=30,
                                 gene_kos=True,
                                 gko_cost=gko_cost,

                                 solution_approach=sd.names.BEST)
# Print solution
print(f"One compressed solution with cost {sols.sd_cost[0]} found and "+\
      f"expanded to {len(sols.gene_sd)} solutions in the uncompressed netork.")
print(f"Example intervention set: {['+' + s if v > 0 else '-' + s for s, v in sols.gene_sd[0].items() if v != 0]}")

import matplotlib.pyplot as plt

# Wild-type plot
datapoints, triang, plot1 = sd.plot_flux_space(ecc_14bdo,
                                               ('BIOMASS_Ecoli_core_w_GAM', 'EX_14bdo_e'),
                                               show=False);
# Plot minimal enforced growth rate
_, _, plot2 = sd.plot_flux_space(ecc_14bdo,
                                 ('BIOMASS_Ecoli_core_w_GAM', 'EX_14bdo_e'),
                                 constraints='BIOMASS_Ecoli_core_w_GAM>=0.5',
                                 show=False);
plot2.set_facecolor('#70AD47')
plot2.set_edgecolor('#70AD47')

# OptKnock design plot
interventions = [[{s: 1.0}, '=', 0.0] for s, v in sols.reaction_sd[0].items() if v < 1]
_, _, plot3 = sd.plot_flux_space(ecc_14bdo,
                                 ('BIOMASS_Ecoli_core_w_GAM', 'EX_14bdo_e'),
                                 # The sign of the glucose exchange reaction is flipped since
                                 # reaction is defined in the direction of secretion.
                                 constraints=interventions,
                                 show=False);
plot3.set_facecolor('#FFC000')
plot3.set_edgecolor('#FFC000')
# adjust axes limits and show plot
plot3.axes.set_xlim(0, 1.05 * max([a[0] for a in datapoints]))
plot3.axes.set_ylim(0, 1.05 * max([a[1] for a in datapoints]))
plt.show()
