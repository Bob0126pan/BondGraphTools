# %% [markdown]
# # Energetic analysis of biochemical systems
# In addition to the flow of mass, energy flows through biochemical networks. Because bond graphs explicily model the flow of power, they naturally enable the analysis of the energetics of biochemical systems. In this notebook we demonstrate how BondGraphTools can be used to examine such energetic quantities. We first create and simulate the Michaelis-Menten model seen in the previous notebook:

# %%
from BondGraphTools import reaction_builder
from BondGraphTools.reaction_builder import Reaction_Network

from BondGraphTools import simulate
from numpy import log, array
import matplotlib.pyplot as plt
from sympy import init_printing, SparseMatrix, Eq, Symbol, lambdify
init_printing()

def initialise_MM_model():
    rn_MM = Reaction_Network(name='Michaelis-Menten enzyme',temperature=310)
    rn_MM.add_reaction('E + S = C', name='R1')
    rn_MM.add_reaction('C = E + P', name='R2')
    rn_MM.add_chemostat('S')
    rn_MM.add_chemostat('P')
    
    model = rn_MM.as_network_model()
    
    (model/"C:E").set_param('k',1)
    (model/"C:C").set_param('k',1)
    (model/"R:R1").set_param('r',1)
    (model/"R:R2").set_param('r',1)
    
    R = reaction_builder.R
    T = 310
    K_S = 1
    K_P = 1
    x_S = 2
    x_P = 1
    (model/"SS:S").set_param('e',R*T*log(K_S*x_S))
    (model/"SS:P").set_param('e',R*T*log(K_P*x_P))
    return rn_MM,model

rn,model = initialise_MM_model()
t,x = simulate(model, timespan=(0.,3.), x0=[1,2])
plt.plot(t,x)
plt.xlabel('Time')
plt.ylabel('Amount')
plt.legend(['E','C'])

# %% [markdown]
# In order to examine the energetics of biochemical systems, we require values (i.e. potentials and mass flows) associated with the bonds of each model rather than the states of the model. These quantities can be extracted using the full equations of the model.
# 
# Below we define a function that outputs the full equations of the model

# %%
def full_equations(model):
    # Load full equations of model
    X, mapping, A, F, G = model.system_model()
    # AX + F(X) = 0
    # G(X) = 0
    AX = A*SparseMatrix(X) + F
    full_model_equations = {}
    for i in range(AX.rows):
        xi = X[i]
        eqn = xi - AX[i,0]
        full_model_equations[str(xi)] = eqn
    return full_model_equations

full_model_equations = full_equations(model)
full_model_equations

# %% [markdown]
# We then search for the relevant bonds connected to the **Re** components and extract the correct mathematical expressions for their fluxes.

# %%
# Define a function that returns the port for a component
def find_port(component,direction):
    if direction in ['f', 'forward']:
        index = 0
    elif direction in ['r', 'reverse']:
        index = 1
    return list(component.ports.keys())[index]

# Returns the mathematical expression for a flux
def reaction_flux_expression(model,Re_comp,direction):
    mapping = model.system_model()[1]
    port = find_port(Re_comp,direction)
    bond_index = mapping[1][port]
    
    full_model_equations = full_equations(model)
    V = full_model_equations[f'f_{bond_index}']
    return V
    
R1 = model/"R:R1"
R2 = model/"R:R2"

V_R1 = reaction_flux_expression(model,R1,'f')
V_R2 = reaction_flux_expression(model,R2,'f')

display(Eq(Symbol("V_1"),V_R1))
display(Eq(Symbol("V_2"),V_R2))

# %% [markdown]
# The flows of the bonds connected to each reaction can be used to calculate their reaction rates. We first extract the corresponding components of the full equation.

# %% [markdown]
# These quantities can be converted into Python functions (taking a state vector as an import) using the `lambdify` function from `sympy`.

# %%
# Returns a function that can be used to compute flux
def flux_function(model,Re_comp):
    V = reaction_flux_expression(model,Re_comp,'f')
    return convert_to_function(V,model)

# Converts a symbolic expression to a function
def convert_to_function(expression,model):
    states = [Symbol(x) for x in model.state_vars.keys()]
    return lambdify(([states]),expression)

V_R1_func = flux_function(model,R1)
V_R2_func = flux_function(model,R2)

# %% [markdown]
# We use the simulation results from the previous section to calculate and plot the reaction rates. The reaction velocities converge to the same value at steady state, as expected.

# %%
V = [[V_R1_func(states), V_R2_func(states)] for states in x]
plt.plot(t,V)
plt.xlabel('Time')
plt.ylabel('Reaction velocity')
plt.legend(['R1','R2'])

# %% [markdown]
# We can also extract the efforts of these bonds, and use them to calculate and plot the affinities of each reaction.

# %%
# Extracts an expression for reaction affinity
def reaction_affinity_expression(model,Re_comp,direction):
    mapping = model.system_model()[1]
    port = find_port(Re_comp,direction)
    bond_index = mapping[1][port]
    
    full_model_equations = full_equations(model)
    A = full_model_equations[f'e_{bond_index}']
    return A

# Extracts the net reaction affinity
def reaction_affinity(model,Re_comp):
    Af = reaction_affinity_expression(model,Re_comp,'f')
    Ar = reaction_affinity_expression(model,Re_comp,'r')
    return Af-Ar

# Extracts a function that can be used to calculate reaction affinity
def reaction_affinity_function(model,Re_comp):
    A = reaction_affinity(model,Re_comp)
    return convert_to_function(A,model)

A_R1 = reaction_affinity_function(model,R1)
A_R2 = reaction_affinity_function(model,R2)

# Calculate the affinities of each and plot
A = [[A_R1(s), A_R2(s)] for s in x]
plt.plot(t,A)
plt.xlabel('Time')
plt.ylabel('Reaction affinity')
plt.legend(['R1','R2'])

# %% [markdown]
# Together with the reaction rates, these give the power consumption of the system.

# %%
power = array(V)*array(A)
total_power = [sum(p) for p in power]
plt.plot(t,total_power)
plt.xlabel('Time')
plt.ylabel('Power')

# %% [markdown]
# It is also possible to extract the chemical potentials of each species using a similar approach. We do so in the code below, and plot how the chemical potentials change with respect to time.

# %%
# Find the component corresponding to a species
def find_species(model,species,metabolites):
    if species in metabolites:
        return model/f'SS:{species}'
    else:
        return model/f'C:{species}'

# Find the port corresponding to a C component
def find_species_port(component,direction):
    return list(component.ports.keys())[0]

# Returns a symbolic expression for the chemical potential of a species
def species_potential(model,species):
    comp = find_species(model,species,metabolites)
    mapping = model.system_model()[1]
    bond_index = mapping[1][find_species_port(comp,0)]
    
    full_model_equations = full_equations(model)
    potential = full_model_equations[f'e_{bond_index}']
    return potential

# Returns a function that can be used to calculate the chemical potential of a species
def species_potential_func(model,species):
    potential = species_potential(model,species)
    return convert_to_function(potential,model)
    
chemical_potentials = {}
metabolites = list(rn._chemostats.keys())

for species in rn.species:
    potential_func = species_potential_func(model,species)
    # Calculate the chemical potential over the simulation
    cp = [potential_func(s) for s in x]
    # Store the results in a dictionary
    chemical_potentials[species] = cp

# Plot the chemical potential for each species
plt.figure()
for species in chemical_potentials.keys():
    plt.plot(t,chemical_potentials[species])
plt.xlabel('Time')
plt.ylabel('Chemical potential')
plt.legend(chemical_potentials.keys())


