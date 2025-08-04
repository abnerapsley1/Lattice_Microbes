# 
# University of Illinois Open Source License
# Copyright 2008-2018 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Michael J. Hallock and Joseph R. Peterson
#
# 

##############################################
### Importing Additional Necessary Modules ###
##############################################
import lm, math, os, h5py
from .LMLogger import *
# Attempt to import tqdm #
try:
	from tqdm import tqdm
# If an error is encountered, define custom tqdm function #
except:
	def tqdm(x,ascii=False):
		return x

##################################
### Define CMESimulation Class ###
##################################
class CMESimulation:
	"""
    This class is the base constructor class for CME simulations.

    Initialization Parameters (Example):
        myObject = jLM.CME.CMESimulation(volume = None, name = "unnamed"):
            volume (float): The volume (in Liters) of the simulation object. If the 
                simulation volume is set to "None", the program assumes that the 
                rate constants are all volume-normalized.
            name (string): The name given to the CMESimulation object.

    Attributes:
        particleMap (dictionary): A dictionary giving order that species were 
            added to the CMESimulation object. The entry is named after the species name
            and the value of the entry is the number representing the order the 
            species was added. This is edited by running the addSpecies method.
        species_id (list): A list of strings giving the names of each chemical 
            species.
        initial_counts (dictionary): A dictionary with the entry being the species
            name and the value being the initial count of the species given in 
            number of particles.
        reactions (list): A list of the reactions in the current model. This list is 
            given as a list of lists in the form ((reactants,products,rate),(reactants,products,rate)...)
        parameters (dictionary): A dictionary containing various simulation parameters 
            such as: seed (for the random number generator), writeInterval (biological time between
            writing results to output file), hookInterval (biological time between communication with
            other simulation objects), maxTime (total biological simulation time), and timestep (amount
            of biological time in a single timestep).
        volume (float): The volume of the reaction simulation (in Liters).
        name (string): The name of the simulation.
        replicates (list): A list of 1 through the number of replicates that have been
            run in simulation for the CMESimulation object.
        filename (string): The name of the file for which simulation results have been saved, if 
            a simulation has been run with the current CMESimulation object.


    Methods: 
        __init__(self, newbase): A dunder method that initializes the instance of the class object.
        defineSpecies(self, species): A method to define the list of chemical species to be used
            in the simulation.
        addParticles(self, species='unknown', count=1): A method to set the initial particle counts
            for a specific species.
        addConcentration(self, species='unknown', conc=0.0): A method to set the initial particle 
            counts by giving a concentration (in Molar) for a specific species.
        addReaction(self, reactant, product, rate): A method to add a reaction with its associated
            reactants, products, and rate constant to the simulation object.
        buildReactionModel(self): A method to convert the current CMESimulation object to an object
            called a ReactionModel from the lm software. This should be run after completely
            specifying the CMESimulation object. The result of this method is used as the input for 
            the simulation methods below. This method does not need to be run independently because
            it is called in the save method below.
        setRandomSeed(self, seed): Set a number for the "seed" entry of the CMESimulation attribute 
            parameter. This determines the random number generator seed for simulation runs.
        setWriteInterval(self, time): Set a number for the "writeInterval" entry of the 
            CMESimulation attribute parameter. This determines the biological time between
            writing results to an output file.
        setHookInterval(self, time): Set a number for the "hookInterval" entry of the 
            CMESimulation attribute parameter. This determines the biological time between 
            communication with other simulation objects.
        setSimulationTime(self, time): Set a number for the "maxTime" entry of the 
            CMESimulation attribute parameter. This determines the total biological time to
            run the current simulation object for each replicate.
        setTimestep(self, time): Set a number for the "timestep" entry of the 
            CMESimulation attribute parameter. This determines the biological time that equates
            to one timestep in the current simulation.
        save(self, filename): A method to save the current CMESimulation object into an HDF5 file
            to be called by the run methods below. 
            Does this method need to be run before any of the run methods below (???)
        run(self, filename, method, replicates = 1, seed = None, cudaDevices = [], checkpointInterval = 0):
            A method to run, in serial, a specified number of simulation replicates with the current 
            CMESimulation object. The solvers available in this method are listed in the method description.
        runMPI(self, filename, method, replicates=1, driver="mpirun", ppe=1, seed=None):
            A method to run, in parallel, a specified number of simulation replicates with the current 
            CMESimulation object. The solvers available in this method are listed in the method description.
        runSolver(self, filename, solver, replicates=1, cudaDevices=None, checkpointInterval=0):
            A method to run, in serial, a specified number of simulation replicates with the current 
            CMESimulation object. This method allows one to specify a user-defined solver.

    """
    
    ### Define Class Initialization Method ###
	def __init__(self,volume = None, name="unnamed"):
        """ 
        Initializer method for the CMESimulation class.
        
        Args:
            self (CMESimulation): Object pointer.
            volume (float): The volume (in Liters) of the simulation object. If the 
                simulation volume is set to "None", the program assumes that the 
                rate constants are all volume-normalized.
            name (string): The name given to the CMESimulation object.
                
        Returns:
            none
        
        """
        
        # Define particleMap attribute as an empty dictionary #
		self.particleMap={}
		# Define species_id attribute as an empty list #
		self.species_id=[]
		# Define initial_counts attribute as an empty dictionary #
		self.initial_counts={}
		# Define reactions attribute as an empty list #
		self.reactions=[]
		# Define parameters attribute as an empty dictionary #
		self.parameters={}
		# Define volume attribute as "volume" argument #
		self.volume = volume
		# Define volume attribute as "name" argument #
		self.name = name
        # Define replicates attribute as an empty list #
		self.replicates = []
		# Set filename attribute to None #
		self.filename = None


    ### Define Class defineSpecies Method ###
	def defineSpecies(self, species):
		""" 
        A method which defines all chemical species in the CMESimulation object.

        This method adds an arbitrary number of chemical species to the CMESimulation.
        It also adds particleMap information for each species and sets the value
        in the initial_counts dictionary for the species to 0.
        
        Args:
            self (CMESimulation): Object pointer.
            species (string or list): A string (when only one species in model) or
                a list (when multiple species in model) that specifies the names of 
                all the chemical species in the CMESimulation object.
                
        Returns:
            none
        
        """

        # Iterate through all entries in "species" argument #
		for s in species:
			# Add the species entry to the species_id attribute #
			self.species_id.append(s)
			# Add a particle map dictionary entry for the current species being added to the CMESimulation object #
			# It seems like the particleMap attribute maps the species number (i.e., how many species before this species + 1) #
			self.particleMap[s]=len(self.species_id)
			# Set the initial_counts dictionary entry of the species being added to 0 #
			self.initial_counts[s]=0




    ### Define Class addParticles Method ###
	def addParticles(self, species='unknown', count=1):
		""" 
        A method which changes initial particle counts for a species.

        This method adds an arbitrary number of chemical species to the CMESimulation.
        It also adds particleMap information for each species and sets the value
        in the initial_counts dictionary for the species to 0.
        
        Args:
            self (CMESimulation): Object pointer.
            species (string): A string (when only one species in model) or
                a list (when multiple species in model) that specifies the names of 
                all the chemical species in the CMESimulation object.
            count (integer): The number of particles to set as the initial_counts 
                attribute for the specified species.
                
        Returns:
            none
        
        """

		try:
			particleNum=self.particleMap[species]
		except KeyError:
			particleNum=1
			LMLogger.warn('In CME.addParticles, couldn\'t find particle of type "%s" in map (is it previously defined with \'CME.defineSpecies(...)\'?).  Shouldn\'t happen.',species)
		self.initial_counts[species]+=count



    ### Define Class addConcentration Method ###
	def addConcentration(self, species='unknown', conc=0.0):
		""" 
        A method which changes initial particle counts for a species.

        This method adds an arbitrary number of chemical species to the CMESimulation.
        It also adds particleMap information for each species and sets the value
        in the initial_counts dictionary for the species to 0.
        
        Args:
            self (CMESimulation): Object pointer.
            species (string): A string (when only one species in model) or
                a list (when multiple species in model) that specifies the names of 
                all the chemical species in the CMESimulation object.
            conc (float): The concentration (in Molar) to add to the initial_counts attribute
                for the specified species. This concentration is then converted to
                a particle count by rounding to the nearest integer using Avogadro's number
                rounded to eight decimal places.
                
        Returns:
            none
        
        """

        # Test whether the volume attribute is set as None #
		if self.volume == None:
			# Print an error using LMLogger that no volume is specified #
			LMLogger.error('In CME.addConcentration, a volume must be specified in CME constructor.')
			raise Exception("No volume specified.")
		# Test whether species is present in particleMap attribute #
		try:
			particleNum=self.particleMap[species]
		# If species is not present, print an error using LMLogger that species is not specified #
		except KeyError:
			particleNum=1
			LMLogger.warn('In CME.addConcentration, couldn\'t find particle of type "%s" in map (is it previously defined with \'CME.defineSpecies(...)\'?).  Shouldn\'t happen.',species)
		# Add species-specific particles to initial_count attribute by converting concentration to particle numbers #
		self.initial_counts[species]= int(round(conc*self.volume*6.02214076e23))



    ### Define Class addReaction Method ###
	def addReaction(self, reactant, product, rate):
		""" 
        A method which adds a reaction to the CMESimulation object.

        This method adds a reaction with its associated reactants, products, and rate
        constant to the CMESimulation object reactions attribute. 
        The reactants and products should be lists and the rate should be positive.
        How do we specify the order of the reaction (zero, first, second)???
        
        Args:
            self (CMESimulation): Object pointer.
            reactant (list): REQUIRED, a list of reactants that participate in the reaction
                being created.
            product (list): REQUIRED, a list of products that are formed in the reaction
                being created.
            rate (float): REQUIRED, the stochastic rte of reaction.
            Should a reaction name be added here as an argument (???)
                
        Returns:
            none
        
        """

        # Test if the rate argument is negative #
		if rate <= 0.0:
			# Return error if rate argument is negative using LMLogger #
			LMLogger.error("In CME.addReaction, the rate must be a positive number")
		# Test if the reaction being added is already in the reactions object attribute #
		if (reactant,product,rate) in self.reactions:
			# Return error if reaction is already in the reactions object attribute #
			LMLogger.warning("Reaction already in model: %s -> %s : %f"%(reactant,product,rate))
			# Exit function #
			return
		# If reaction is not already in reactions object attribute, add the reaction as a list #
		# containing the reactant list, the product list, and the rate constant. #
		self.reactions.append((reactant,product,rate))




    ### Define Class buildReactionModel Method ###
	def buildReactionModel(self):
		""" 
        A method which builds a lattice microbes ReactionModel from a CMESimulation object.

        This method adds a reaction with its associated reactants, products, and rate
        constant to the CMESimulation object reactions attribute. 
        The reactants and products should be lists and the rate should be positive.
        This method should be used AFTER defining the spcies, initial abundances, 
        reactions, parameters, and the reactants and products for each reaction in the 
        CMESimulation object (i.e., defineSpecies, addParticles, and addReaction methods).
        
        Args:
            self (CMESimulation): Object pointer.
                
        Returns (ReactionModel from lm):
            A ReactionModel object as defined in lm.
        
        """

        # Initialize rm as a ReactionModel object from the lm modules #
		rm=lm.ReactionModel()

        # Define the number of reactions in input CMESimulation object #
		numReactions=len(self.reactions)
		# Add number of reactions in CMESimulation object to info file using LMLogger #
		LMLogger.info("number of reactions = %d", numReactions)
		# Set the number of reactions of rm to be equa to the CMESimulation reaction number #
		rm.set_number_reactions(numReactions)

        # Define the number of species in input CMESimulation object #
		numSpecies=len(self.species_id)
		# Add number of species in CMESimulation object to info file using LMLogger #
		LMLogger.info("number of species = %d", numSpecies)
		# Set the number of species of rm to be equal to the CMESimulation species number #
		rm.set_number_species(numSpecies)

        # Iterate through all species in CMESimulation object #
		for s in self.species_id:
			# Add information on initial particle counts for species to info file using LMLogger #
			LMLogger.debug("\t set initial counts of %s (id %d) to %d", s, self.particleMap[s], self.initial_counts[s])
			# Set the initial count of particles for species in rm object to be equal to species initial_counts attribute in CMESimulation object #
			rm.add_initial_species_count(self.initial_counts[s])

		# Initialize reaction lists/matrices #
		# Define rxtypes list as a list of 0s with numReactions length #
		# This will be a list of reaction types #
		# Zero-order: 0; First-order: 1; Second-order: 2; Second-order self reaction: 3 #
		rxtypes = [0] * numReactions
		# Define rxconst list as a list of 0s with numReactions length #
		# This will be a list of reaction constants for each reaction #
		rxconst = [0] * numReactions
		# Define stoich list of length numSpecies * numReactions #
		# This will store stoichiometry information (-1 for reactant, 1 for product, -2 for second-order self-reaction reactant) 
		# for each species in each reaction #
		# Can be thought of as a pseudo-matrix or list of lists, even though it is just a long list #
		stoich  = [0] * numSpecies * numReactions
		# Define depmat matrix as a list of numReactions length of lists of 0s with numSpecies length #
		# This will store dependancy information (entry will be 1 if species participates in reaction) #
		# for each species in each reaction #
		# Can be thought of as a pseudo-matrix or list of lists, even though it is just a long list #
		depmat  = [0] * numSpecies * numReactions

        # Initialize rnum (reaction number) variable to iterate through all reactions #
		rnum=0
		# Iterate through all reactions attributes in CMESimulation object #
		for rx in self.reactions:
			# Define reactant list #
			reactant=rx[0]
			# Define product list #
			product=rx[1]
			# Define rate constant #
			rate=rx[2]

			# Change rates if CMESimulation object has a specified volume #
			# Test if volume attribute is not set as None #
			if self.volume != None:
				# Test if reactant is a tuple (multiple values) #
				if isinstance(reactant, tuple):
					# Test if the reactants are different (second-order reaction) #
					if reactant[0] != reactant[1]:
						# Divide the rate constant by Avogadro's number and volume attribute #
						rate /= (6.022e23*float(self.volume))
					# Test if the reactants are the same (second-order self-reaction) #
					elif reactant[0] == reactant[1]:
						# Divide the rate constant by Avogadro's number and volume attribute #
						rate /= (6.022e23*float(self.volume))
				# If there are one or no reactants for the reaction... #
				else:
					# Test if there are no reactants (zero-order reaction) #
					if reactant == '':
						# Multiply the rate constant by Avogadro's number and volume attribute #
						rate *= 6.022e23*float(self.volume)
					# Test if there is one reactant (first-order reaction) #
					# do nothing...

            # Write reaction information to debug file using LMLogger #
			LMLogger.debug("\t rx %d is %s -> %s at rate %g", rnum, reactant, product, rate)
			
			# Add the rate constant for the specified reaction to the rxconst list #
			rxconst[rnum]=rate

			# Add the reaction order information to the rxtypes list #
			# Test if reactants for reaction are a tuple (i.e., second-order) #
			if(isinstance(reactant, tuple)):
				# Test if reactants are not equal (second-order reaction)#
				if reactant[0] != reactant[1]:
					# Set the rxtypes equal to 2 for specified reaction #
					rxtypes[rnum]=2
					# Iterate through both reactants #
					for sr in reactant:
						# Update depmat entry for the specific reactant in the reaction #
						depmat[rnum + numReactions*(self.particleMap[sr]-1)] = 1
						# Update stoich entry for the specific reactant in the reaction #
						stoich[rnum + numReactions*(self.particleMap[sr]-1)] -= 1
				# Test if reactants are equal (second-order self reaction) #
				elif reactant[0] == reactant[1]:
				    # Set the rxtypes equal to 3 for specified reaction #
					rxtypes[rnum]=3
					# Update depmat entry for the specific reactant in the reaction #
					depmat[rnum + numReactions*(self.particleMap[reactant[0]]-1)] = 1
					# Update stoich entry for the specific reactant in the reaction #
					stoich[rnum + numReactions*(self.particleMap[reactant[0]]-1)] -= 2
			# All reactions that are either zero- or first-order #
			else:
				# Test if there are no reactants (zero-order) #
				if reactant == '':
					# Set the rxtypes equal to 0 for specified reaction #
					rxtypes[rnum] = 0
					# No need to update stoich entry for the specific reactant in the reaction (0) #
				# Test if there is one reactant (first-order) #
				else:
				    # Set the rxtypes equal to 1 for specified reaction #
					rxtypes[rnum]=1
					# Update depmat entry for the specific reactant in the reaction #
					depmat[rnum + numReactions*(self.particleMap[reactant]-1)] = 1
					# Update stoich entry for the specific reactant in the reaction #
					stoich[rnum + numReactions*(self.particleMap[reactant]-1)] -= 1
            
            # Determine if there are multiple products for specified reaction # 
			if(isinstance(product, tuple)):
				# Iterate through products #
				for sp in product:
					# Update stoich entry for the specific product in the reaction #
					stoich[rnum + numReactions*(self.particleMap[sp]-1)] += 1
			# Determine if there is one or zero products for specified reaction #
			else:
				# Test if there is one product in the reaction #
				if product != '':
					# Update stoich entry for the specific product in the reaction #
					stoich[rnum + numReactions*(self.particleMap[product]-1)] += 1
            
            # Increment rnum to move on to next reaction #
			rnum += 1


        # Print reaction type to debug file using LMLogger #
		LMLogger.debug("rxtypes: %s", rxtypes)
		# Print reaction constant to debug file using LMLogger #
		LMLogger.debug("rxconst: %s", rxconst)

		# Iterate through all reactions in CMESimulation object #
		for r in range(numReactions):
			# Add a reaction to the lm ReactionModel object #
			rm.add_reaction()
			# Add reaction type to the lm ReactionModel object #
			rm.mutable_reaction(r).set_type(rxtypes[r])
			# Add rate constant to the lm ReactionModel object #
			rm.mutable_reaction(r).add_rate_constant(rxconst[r])

        # Print reaction stoichiometry list in debug file using LMLogger #
		LMLogger.debug("stoich:  %s", stoich)
		# Iterate through all entries of stoich #
		for x in range(len(stoich)):
			# Add stoichiometry information to lm ReactionModel object #
			rm.add_stoichiometric_matrix(stoich[x])

        # Print reaction dependency list in debug file using LMLogger #
		LMLogger.debug("depmat:  %s", depmat)
		# Iterate through all entries of depmat #
		for x in range(len(depmat)):
			# Add dependency information to lm ReactionModel object #
			rm.add_dependency_matrix(depmat[x])

        # Return the new ReactionModel object from lm #
		return rm




    ### Define Class setRandomSeed Method ###
	def setRandomSeed(self, seed):
		""" 
        A method which sets a random seed for the CMESimulation object.
        
        Args:
            self (CMESimulation): Object pointer.
            seed (int or string): REQUIRED, a selected number for the initialization of the 
                random number generator.
                
        Returns:
            none
        
        """

        # Add seed input argument to the parameters attribute dictionary and name it "seed" #
		self.parameters['seed']=str(seed)





    ### Define Class setWriteInterval Method ###
	def setWriteInterval(self, time):
		""" 
        A method which sets a write-to-disk interval for the simulation.

        The interval specified will determine the (biological) time that must 
        pass before the simulation writes the state of the cell to an output file.
        
        Args:
            self (CMESimulation): Object pointer.
            time (int or string): REQUIRED, a selected number for the time between
                writing events during the simulation.
                
        Returns:
            none
        
        """

        # Add time input argument to the parameters attribute dictionary and name it "writeInterval" #
		self.parameters['writeInterval']=str(time)




    ### Define Class setHookInterval Method ###
	def setHookInterval(self, time):
		""" 
        A method which sets a whook interval for the simulation.

        The interval specified will determine the (biological) time that must 
        pass before the simulation communicates the state of the cell to 
        other simulation objects. This is used for integrating various types of
        simulations into one model (i.e., combining CME and ODE).
        
        Args:
            self (CMESimulation): Object pointer.
            time (int or string): REQUIRED, a selected number for the time between
                hook events during the simulation.
                
        Returns:
            none
        
        """

        # Add time input argument to the parameters attribute dictionary and name it "hookInterval" #
		self.parameters['hookInterval']=str(time)




    ### Define Class setSimulationTime Method ###
	def setSimulationTime(self, time):
		""" 
        A method which sets a total simulation time.

        The interval specified will determine the total (biological) simulation time.
        
        Args:
            self (CMESimulation): Object pointer.
            time (int or string): REQUIRED, a selected number for the total duration
                in biological time for the simulation.
                
        Returns:
            none
        
        """

        # Add time input argument to the parameters attribute dictionary and name it "maxTime" #
		self.parameters['maxTime']=str(time)




    ### Define Class setTimestep Method ###
	def setTimestep(self, time):
		""" 
        A method which sets a simulation timestep.

        The interval specified will determine the amount of (biological) time that must pass
        during the simulation to be considered one time-step.
        
        Args:
            self (CMESimulation): Object pointer.
            time (int or string): REQUIRED, a selected number for the total duration
                in biological time for the simulation timestep.
                
        Returns:
            none
        
        """

        # Add time input argument to the parameters attribute dictionary and name it "timestep" #
		self.parameters['timestep']=str(time)




    ### Define Class save Method ###
	def save(self, filename):
		""" 
        A method which saves the current CMESimulation object configuration to an HDF5 file.

        This method can be used to save the current model configuration as a lattice for either bookkeeping
        or use in subsequent runs.
        Does this method need to be run prior to running the run/solve methods (???)
        
        Args:
            self (CMESimulation): Object pointer.
            filename (string): REQUIRED, a name for the file that will store 
                the current .
                
        Returns:
            none
        
        """

        # Use lm SimulationFile.create() method to create a file named after the filename input argument #
		lm.SimulationFile.create(filename)
		# Define f as a lm SimulationFile object initialized with the filename input argument #
		f=lm.SimulationFile(filename)

        # Define rm as the lm ReactionModel output of the current CMESimulation object #
		rm=self.buildReactionModel()
		# Use the ReactionModel rm as input for the lm SimulationFile.setReactionModel() method #
		# This will push the lm ReactionModel object into the newly created HDF5 file #
		f.setReactionModel(rm)

        # Iterate through all items in the CMESimulation parameters attribute (a dictionary) #
		for key in self.parameters:
			# Print parameter being set and its value to debug file using LMLogger #
			LMLogger.debug("Setting parameter %s = %s", key, self.parameters[key])
			# Use lm SimulationFile.setParameter() method to save parameter from CMESimulation object to #
			# the lm ReactionModel in the HDF5 file #
			f.setParameter(key, self.parameters[key])

		# Print all species names in CMESimulation object to debug file using LMLogger #
		LMLogger.debug("SpeciesNames: %s", ",".join(self.species_id))

        # Use lm SimulationFile.close() method to close the new HDF5 file #
		f.close()
    
        # Define the simFile as an h5py.File object which contains the newly created HDF5 file #
        # Enable read/write capabilities for this new HDF5 file #
		simFile = h5py.File(filename,'r+')
    
        # Define dt using the h5py.special_dtype() function #
        # This function enables the creation of NumPy dtype objects that do not have #
        # native equivalents in NumPy. Here, we use the vlen (variable-length) data type #
        # where individual elements within the dataset can have varying lengths, such #
        # as lists of strings. We define the base type in our list to be a string #
		dt = h5py.special_dtype(vlen=str)
		# Iterate through the CMESimulation object attribute species_id and perform #
		# the encode() function on each species name. The encode function will encode the #
		# species names into ASCII format and ignore any characters that are not included #
		# in the ASCII format (e.g., Unicode text). #
		# The result of this iteration is a list of species names stored in the object #
		# specNamesListAscii #
		specNamesListAscii = [n.encode("ascii", "ignore") for n in self.species_id]
        # Save the specNamesListAscii object in the new HDF5 file with simulation data #
        # This species list will be saved in the "Parameters/SpeciesNames" location #
		simFile.create_dataset('Parameters/SpeciesNames', (len(specNamesListAscii),1), dt, specNamesListAscii)
        # Force HDF5 data from buffer to disk using the h5py.File.flush() method #
		simFile.flush() 
		# Close the new HDF5 file using the h5py.File.close() method #
		simFile.close() 




    ### Define Class run Method ### 
	def run(self, filename, method, replicates = 1, seed = None, cudaDevices = [], checkpointInterval = 0):
		""" 
        A method which runs the current CMESimulation object simulation.

        This method is used to perform the simulation run on the current CMESimulation object.
        This method results in a new file being created (filename argument) that stores the
        data from the simulation run(s). The stochastic solver method can be specified
        along with the number of simulation replicates, the seed for the random number generator,
        which GPUs are to be used, and the checkpointInterval.
        
        Args:
            self (CMESimulation): Object pointer.
            filename (string): REQUIRED, a name for the file that will store 
                the current simulation run output.
            method (string): REQUIRED, the method of stochastic solver to be used
                for the present simulation runs. CME options include:
                    lm::cme::FluctuatingNRSolver
                    lm::cme::GillespieDSolver
                    lm::cme::HillSwitch
                    lm::cme::LacHillSwitch
                    lm::cme::NextReactionSolver
                    lm::cme::SelfRegulatingGeneSwitch
                    lm::cme::TwoStateExpression
                    lm::cme::TwoStateHillSwitch
                    lm::cme::TwoStateHillLoopSwitch
            replicates (integer): The number of simulation replicates to run.
            seed (string or float): The seed for the random number generator. This could 
                have been set prior to running the simulations using the setRandomSeed
                method from the CMESimulation class. If a value is specified, it will
                take precedence over any previously defined seed value.
            cudeDevices (list): An list of integers that represent the GPU-CUDA 
                devices to be used in the simulation (???).
            checkpointInterval (integer): An interval in time which, if passed, will
                run the lm software lm::main::CheckpointSignaler(). If set at 0,
                no checkpointIntervals will be run. 

        Returns:
            none
        
        """

        # Test whether input argument seed is not set as None #
		if seed is not None:
			# Define f as an lm SimulationFile object instantiated with #
			# the filename input argument (new HDF5 file to be generated) #
			f = lm.SimulationFile(filename)
			# Run lm SimulationFile.setParameter method to save input argument seed #
			# as a parameter in the output file named "seed" #
			f.setParameter("seed", str(seed))
			# Close new output file using lm SimulationFile.close method #
			f.close()

        # Iterate from 1 to the number of replicates defined at input #
		for r in tqdm(range(1, replicates + 1), ascii = True):
			# Print message that you are running custom replicate simulations to debug using LMLogger #
		    LMLogger.debug("Running replicate %d% with Solver: %s"%(r,solver))
			# Run the specified replicate using lm class run.Simulation method #
			# Input for the runSimulation method arguments comes from the 
			# previously specified arguments for the run method #
			# r is the replicate number being run #
			lm.runSimulation(filename, r, method, cudaDevices, checkpointInterval)
			# Update internal state with replicates that have been run by #
			# appending the most recently run replicate number to the CMESimulation # 
			# object attribute replicates #
			self.replicates.append(r)
			
		# Update the CMESimulation object attribute filename  with input filename argument #
		# This is the file that has had results printed to it #
		self.filename = filename


	

    ### Define Class runMPI Method ###
	def runMPI(self, filename, method, replicates=1, driver="mpirun", ppe=1, seed=None):
		""" 
        A method which runs the current CMESimulation object simulation in parallel.

        This method is used to perform the simulation run on the current CMESimulation object
        using a parallel processing method.
        This method results in a new file being created (filename argument) that stores the
        data from the simulation run(s). The stochastic solver method can be specified
        along with the number of simulation replicates, the software used to run simulations in
        parallel, the number of processing elements (i.e., nodes), and the the seed for the 
        random number generator.
        
        Args:
            self (CMESimulation): Object pointer.
            filename (string): REQUIRED, a name for the file that will store 
                the current simulation run output.
            method (string): REQUIRED, the method of stochastic solver to be used
                for the present simulation runs. CME options include:
                    lm::cme::FluctuatingNRSolver
                    lm::cme::GillespieDSolver
                    lm::cme::HillSwitch
                    lm::cme::LacHillSwitch
                    lm::cme::NextReactionSolver
                    lm::cme::SelfRegulatingGeneSwitch
                    lm::cme::TwoStateExpression
                    lm::cme::TwoStateHillSwitch
                    lm::cme::TwoStateHillLoopSwitch
            replicates (integer): The number of simulation replicates to run.
            driver (string): The name of the software being used to run simulations in parallel.
                Options include (but are not limited to): mpirun, aprun, and ibrun.
            ppe (integer): The number of parallel processing elements to use (i.e., the number
                of nodes for the simulation).
            seed (string or float): The seed for the random number generator. This could 
                have been set prior to running the simulations using the setRandomSeed
                method from the CMESimulation class. If a value is specified, it will
                take precedence over any previously defined seed value.

        Returns:
            none
        
        """

		# Test whether input argument seed is not set as None #
		if seed is not None:
			# Define f as an lm SimulationFile object instantiated with #
			# the filename input argument (new HDF5 file to be generated) #
			f = lm.SimulationFile(filename)
			# Run lm SimulationFile.setParameter method to save input argument seed #
			# as a parameter in the output file named "seed" #
			f.setParameter("seed", str(seed))
			# Close new output file using lm SimulationFile.close method #
			f.close()

		# Define repStr for storing the number of replicates to simulate #
		repStr = "1"
		# Test if replicates argument input is not equal to 1 #
		if replicates != 1:
			# Redefine repStr to be "1-#reps" #
			repStr = "1-%d"%replicates
		# Define a string called cmdStr for storing the command line executable code #
		# All variables are taken from input arguments or repStr #
		# driver software is used to run the Lattice Microbes from the command line #
		cmdStr = '%s -np %d  lm -ws -r %s -sl "%s"  -f %s'%(driver, ppe, repStr, method, filename)
		# Print which command was run to the debug file using LMLogger #
		LMLogger.debug("Running MPI LM job with command:\n\t%s"%(cmdStr))

		# Run the cmdStr command and save the return value as status (should be 0). #
		status = os.system(cmdStr)
		
		# Test if there was a non-zero return status from the previous command #
		if status != 0:
			# Print an error associated with the cmdStr to error file using LMLogger #
			LMLogger.error("Failed running MPI (status: %d) job running command:\n\t%s\n"%(status, cmdStr))

		# Iterate from 1 to the number of replicates defined at input #
		for r in range(1,replicates):
			# Update internal state with replicates that have been run by #
			# appending the most recently run replicate number to the CMESimulation # 
			# object attribute replicates #
			self.replicates.append(r)

		# Update the CMESimulation object attribute filename  with input filename argument #
		# This is the file that has had results printed to it #
		self.filename = filename




    ### Define Class runSolver Method ###
	def runSolver(self, filename, solver, replicates=1, cudaDevices=None, checkpointInterval=0):
		""" 
        A method which runs the current CMESimulation object simulation in serial with a custom solver.

        This method is used to perform the simulation run on the current CMESimulation object
        in serial using a custom solver.
        This method results in a new file being created (filename argument) that stores the
        data from the simulation run(s). The stochastic solver method can be specified
        along with the number of simulation replicates, which GPUs are to be used, and the 
        checkpointInterval. Should seed be an argument here also (???)
        
        Args:
            self (CMESimulation): Object pointer.
            filename (string): REQUIRED, a name for the file that will store 
                the current simulation run output.
            solver (string): REQUIRED, the method of stochastic solver to be used
                for the present simulation runs. this is a custom solver input:
                    lm::cme::FluctuatingNRSolver
                    lm::cme::GillespieDSolver
                    lm::cme::HillSwitch
                    lm::cme::LacHillSwitch
                    lm::cme::NextReactionSolver
                    lm::cme::SelfRegulatingGeneSwitch
                    lm::cme::TwoStateExpression
                    lm::cme::TwoStateHillSwitch
                    lm::cme::TwoStateHillLoopSwitch
            replicates (integer): The number of simulation replicates to run.
            cudeDevices (list): An list of integers that represent the GPU-CUDA 
                devices to be used in the simulation. ???
            checkpointInterval (integer): An interval in time which, if passed, will
                run the lm software lm::main::CheckpointSignaler(). If set at 0,
                no checkpointIntervals will be run. 

        Returns:
            none
        
        """

        # Test if cudeDevices input argument is set as None #
		if cudaDevices is None:
			# Define cudaDevices as a list containing a single 0 #
			# Should this list be empty, or does it matter (???) #
			cudaDevices = [0]
	    # Print message that you are running custom solver simulations to debug using LMLogger #
		LMLogger.debug("Running Custom CMESolver: %s"%(solver))
		# Iterate from 1 to the number of replicates defined at input #
		for r in tqdm(range(1, replicates+1),ascii=True):
			# Print message that you are running a specific replicate to debug using LMLogger #
			LMLogger.debug("  Running replicate: %d"%(r))
			# Run the specified replicate using lm class run.Simulation method #
			# Input for the runSimulation method arguments comes from the 
			# previously specified arguments for the runSolver method #
			# r is the replicate number being run #
			lm.runSolver(filename, r, solver, cudaDevices, checkpointInterval)
			# Update internal state with replicates that have been run by #
			# appending the most recently run replicate number to the CMESimulation # 
			# object attribute replicates #
			self.replicates.append(r)

		# Update the CMESimulation object attribute filename  with input filename argument #
		# This is the file that has had results printed to it #
		self.filename = filename
