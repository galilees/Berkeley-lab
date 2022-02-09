from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import make_sub_header
#%%
import matplotlib as plt
plt.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
from matplotlib import pyplot as plt
import numpy as np

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
mode = queue *one_cpu test
  .type = choice
  .help = Run on the queue or one by one (one processor)
models
  .multiple=true
  {
  pdb_code = None
  .type = str
  }
'''

# =============================================================================

class AnalyseGol(ProgramTemplate):
  description = '''
Script to do xxxxxx
'''
  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  #-----------------------------------------------------------------------------

  def validate(self):
    '''
    Validate inputs
    '''
    make_sub_header('Validating inputs', out=self.logger)
    print('Running in mode: ', self.params.mode, file=self.logger)
    self.data_manager.has_models(raise_sorry=True)

  #-----------------------------------------------------------------------------

  def run(self):
    '''
    Code below will be executed
    '''
    print('here we go', file=self.logger)
    #
    if self.params.mode == 'test':
      self.perform_tests()
      return
    #
    self.model = self.data_manager.get_model()
    self.get_gol_selection()
    self.count_nearby_GOL()
    self.plot_counts()

  #-----------------------------------------------------------------------------

  def get_gol_selection(self):
    '''
    Prints the selection string and iselection for each GOL
    '''
    make_sub_header('Getting selection for GOL', out=self.logger)
    self.gol_selection_dict ={}

    hierarchy = self.model.get_hierarchy()
    for m in hierarchy.models():            # Get hierarchy object
      for chain in m.chains():              # loop over chain, residue group, and atom group 
        for rg in chain.residue_groups():
          for ag in rg.atom_groups():
            
            if (ag.resname == "GOL"):      
            
              iselection = ag.atoms().extract_i_seq()
              sel_str_ligand = " ".join(['chain', chain.id, 'and resname', ag.resname, 'and resseq', rg.resseq])
             
              self.gol_selection_dict[sel_str_ligand] = iselection

    print(self.gol_selection_dict)
              
    
   
    
#-----------------------------------------------------------------------------
  # def nearby_residues_selection_list(self):


  #   b_max = 100
  
    
              
  #             atoms = ag.atoms()
  #             ave = 0 
            
  #             for a in atoms:
  #               ave = ave + a.b 
  #             ave = ave/len(atoms) 
      
            
  #             if ave > b_max:
  #               continue
                
  #             selection_list.append(sel_str_ligand)
            
  
         

#-----------------------------------------------------------------------------
  def count_nearby_GOL(self):
    make_sub_header('Getting residue counts ditionary', out=self.logger)
    resname_dict_list = []


    for sel_str in self.gol_selection_dict.keys():

      near = 'residues_within(5,%s)'%sel_str
      selection_bool1 = self.model.selection(near)
      m1 = self.model.select(selection_bool1)        
      ph1 = m1.get_hierarchy()
      dict_neighbors = ph1.overall_counts().resnames
      
      resname_dict_list.append(dict_neighbors)
    
    res_list = ["Other","HOH","GOL","LYS", "ASN", "TRP", "ASP","GLU","GLN","ALA","SER","HIS","THR","TYR","PHE","ARG","VAL","GLY","CYS","PRO","ILE","MET"]
    self.res_dict = {}
    for aa in res_list:
        self.res_dict[aa] = 0

    for self.resname_dict in resname_dict_list:
      for k, i in self.resname_dict.items():
        for aa in res_list:
          if k == aa:
            self.res_dict[aa] += i
        if k not in res_list:
          self.res_dict["Other"] += i



    print(self.res_dict)
   
   

#-----------------------------------------------------------------------------
  def plot_counts(self):
   
    plt.bar(list(self.res_dict.keys()), self.res_dict.values(), color='r')
    plt.title('Nearby Amono Acids Count', y=1.16, fontsize=14)

    plt.ylabel('Number Amino Acids')
    plt.xlabel('Name of Amino Acid')

    spacing = 0.500
   
    plt.close()
    return
  plt.show(block=False)
#-----------------------------------------------------------------------------

def perform_tests(self):
  '''
  Here we could do tests.
  '''
  make_sub_header('Performing tests', out=self.logger)

# =============================================================================

if __name__ == '__main__':
  #
  from iotbx.cli_parser import run_program
  run_program(program_class=AnalyseGol, args=["/Users/GalileeS/Desktop/test/1bg4.pdb"] )