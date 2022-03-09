from __future__ import division, print_function
from xml.etree.ElementTree import iselement
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
from libtbx.str_utils import make_sub_header

import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from libtbx.utils import Sorry, null_out

import matplotlib as plt
import json, os, sys
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
    self.success = True
    #
    print('Using file %s' % self.data_manager.get_default_model_name())
    #
    self.initialize()
    #
    if self.params.mode == 'test':
      self.perform_tests()
      return
    self.json_data = {}
    #
    self.model = self.data_manager.get_model()
    try:
      self.get_selection("GOL")
    except Exception as e:
      self.success   = False
      print('failed to get GOL selection.\n' , file=self.logger)
       
      self.validate_gol()
    self.res_nearby_count()
    self.save_json
    self.get_hbonds()
   
    # self.count_nearby() # overall count of residues near selection for entire protein
    # self.plot_counts()
    # #self.ave_resdict_aa_dict()
    # self.max_min_res()
    # #self.aa_dict()


  #-----------------------------------------------------------------------------

  def initialize(self):
    '''
    inistialze data structures.
    '''
    self.selection_dict = {}
    self.nearby_residue_dict = {}
    self.res_dict = {}

#-----------------------------------------------------------------------------
  def save_json(self):
      self.json_data['success'] = self.success
      basename, extension = os.path.splitext(os.path.basename(self.data_manager.get_default_model_name()))
      json_filename = basename + '_dataGOL.json'
    
      with open(json_filename, 'w') as fp:
        json.dump(self.json_data, fp, sort_keys=True, indent=4)


  #-----------------------------------------------------------------------------

  def get_hbonds(self):
    '''
    Searches within a radius of 5 angstroms for hydrogen bonds to glycerol and returns a list of tuples containing the isequences of hbond partners 
    '''
    make_sub_header('H-bonds', out=self.logger)
    iselection_list = []
    hbonds_List=[]
    iselection_dict = {}
    for sel_str in self.selection_dict.keys():
      print('Now looking at ', sel_str)
      near_res = 'residues_within(5,%s)'%sel_str
      selection_bool2 = self.model.selection(near_res)
      m2 = self.model.select(selection_bool2)
      m2.set_log(log = null_out())
      m2.process(make_restraints=True)

      pnps = pnp.manager(model = m2)
    
      hbonds = pnps.get_hbonds()
      hbonds.show(log=sys.stdout) 
      print("tuple keys to hbonds table: ", hbonds._hbonds_dict.keys()) #Show tuple index to table
      hierarchy = m2.get_hierarchy()

      for m in hierarchy.models():            # Get hierarchy object
        for chain in m.chains():              # loop over chain, residue group, and atom group
          for rg in chain.residue_groups():
            for ag in rg.atom_groups():
              if (ag.resname == "GOL"):
                iselection = ag.atoms().extract_i_seq()
                iselection_dict[sel_str] = iselection
      for k,v in iselection_dict.items():
        iselection_list.append(list(v))
        iselection_list_flat = list(np.concatenate(iselection_list).flat)
      
      for iseq_tuple in list(hbonds._hbonds_dict.keys()):
        for iseq in iseq_tuple:
          if iseq in iselection_list_flat:
            if iseq_tuple not in hbonds_List:
              hbonds_List.append(iseq_tuple)
    print("list of tuples containing isequences of glycerols with hbonds: ", hbonds_List)
    
 #-----------------------------------------------------------------------------  

  def get_selection(self, residue):
    '''
    residue:str
    takes in the name of a residue from the pdb and 
    Prints a dictionary with the selection string as the key and iselection for the residue as the value.
    '''
    make_sub_header('Getting selection for residue', out=self.logger)

    print(residue)
    hierarchy = self.model.get_hierarchy()
    for m in hierarchy.models():            # Get hierarchy object
      for chain in m.chains():              # loop over chain, residue group, and atom group 
        for rg in chain.residue_groups():
          for ag in rg.atom_groups():
            
            if (ag.resname == residue):      
            
              iselection = ag.atoms().extract_i_seq()
              sel_str = " ".join(['chain', chain.id, 'and resname', ag.resname, 'and resseq', rg.resseq])
             
              self.selection_dict[sel_str] = iselection

              
      self.json_data['selection_strings'] = list(self.selection_dict.keys())
      self.save_json()
    
    for sel_str in self.selection_dict.keys():
      print(sel_str)
    
              
#----------------------------------------------------------------------------
  def validate_gol(self, b_max = 100,occ_min = .2):
    '''
    b_max: int
    occ_min: int
    removes strutures from the selection dictionary that do no meet criteria. The b factor is out of range if it exceeds b_max and the occupancy is out of range if it is less than occ_min. 
    Also checks that all occupancies within the molecule are the same.
    '''
    make_sub_header('curate gol selection', out=self.logger)
    
    bad_selection=[]
   
    for sel_str in self.selection_dict.keys():
      selection_bool1 = self.model.selection(sel_str)
      m1 = self.model.select(selection_bool1)        
      ph1 = m1.get_hierarchy()

      occ_list = list(m1.get_atoms().extract_occ())
      mmm = m1.get_atoms().extract_occ().min_max_mean()
      if occ_list.count(mmm.min != mmm.max):
        bad_selection.append(sel_str)
      

      mmm = m1.get_atoms().extract_occ().min_max_mean()
      if  mmm.mean < occ_min:
       bad_selection.append(sel_str)

      #b_max = 100.0
      mmm = m1.get_atoms().extract_b().min_max_mean()
      if  mmm.mean > b_max:
       bad_selection.append(sel_str)
      
    for s in bad_selection:
      if s in self.selection_dict:
       self.selection_dict.pop(s)
    
    if bad_selection:
      print("selections removed: ")
      for sel_str in bad_selection:
        print(sel_str)
    else:
      print("All Gol passed the curation", file=self.logger)
    
#----------------------------------------------------------------------------  

  def res_nearby_count(self):
    '''
   counts the residues nearby each selection. Returns a dictionary with the name of the residue as the key and its count as the value. 
    '''
    make_sub_header('count of residues nearby each selection', out=self.logger)
    
    
    for sel_str in self.selection_dict.keys():
     
      near_res = 'residues_within(5,%s)'%sel_str
      selection_bool2 = self.model.selection(near_res)
      m = self.model.select(selection_bool2) 
      ph = m.get_hierarchy()
      res_nearby = ph.overall_counts().resnames
      print(res_nearby, file=self.logger)
      self.nearby_residue_dict[sel_str] = res_nearby
      #self.json_data['nearby_res'] = res_nearby
      #self.save_json()

    self.json_data['nearby_res'] = self.nearby_residue_dict
    self.save_json()
  
    #return res_nearby

  #----------------------------------------------------------------------------  

  def count_nearby(self):
    '''
    counts the number of residues nearby a selected ligand and returns a dictionary of the overall counts of residues nearby the selection for the entire protein. 
    '''  
    make_sub_header('Getting residue counts ditionary', out=self.logger)
    resname_dict_list = []


    for sel_str in self.selection_dict.keys():

      near = 'residues_within(5,%s)'%sel_str
      selection_bool1 = self.model.selection(near)
      m1 = self.model.select(selection_bool1)        
      ph1 = m1.get_hierarchy()
      dict_neighbors = ph1.overall_counts().resnames
      
      resname_dict_list.append(dict_neighbors)
    
    res_list = ["Other","HOH","GOL","LYS", "ASN", "TRP", "ASP","GLU","GLN","ALA","SER","HIS","THR","TYR","PHE","ARG","VAL","GLY","CYS","PRO","ILE","MET"]

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
    '''
    plots the counts of residues nearby a selection
    '''  
    
    plt.bar(list(self.res_dict.keys()), self.res_dict.values(), color='r')
    plt.title('Nearby Amono Acids Count', y=1.16, fontsize=14)

    plt.ylabel('Number Amino Acids')
    plt.xlabel('Name of Amino Acid')

    spacing = 0.500
    
    plt.savefig('bla.png')
    plt.close()
    


  #-----------------------------------------------------------------------------
  # def aa_dict(self):
  #   '''
  #   counts the number of amino acids nearby GOL
  #   ''' libtbx.python  /Berkeley-lab/run.py 1bg4.pdb
  #   make_sub_header('Getting amino acid counts ditionary', out=self.logger)  
  #   self.aa_dict = {k: v for k, v in self.res_dict.items() }
  #   self.aa_dict.pop('HOH')
  #   self.aa_dict.pop('GOL')li
  #   self.aa_dict.pop('Other')
  #   print(self.aa_dict)

  #-----------------------------------------------------------------------------
  def ave_resdict_aa_dict(self):
    '''
    returns dictionaries of the average number of residues and the number of amino acids nearby a selected ligand.
    '''
    make_sub_header('Getting average residue counts ditionary and the number of amino acids nearby selected ligand', out=self.logger)  
    ave_resdict = {k: v/len(self.selection_dict) for k, v in self.resname_dict.items()}

    self.aa_dict = {k: v for k, v in self.res_dict.items() }
    self.aa_dict.pop('HOH')
    self.aa_dict.pop('GOL')
    self.aa_dict.pop('Other')
    
    print("Average residues nearby", ave_resdict," Amino acids nearby: ", self.aa_dict) 
  #-----------------------------------------------------------------------------
  def max_min_res(self):
    '''
    maximum and minimum of residues in the ooverall count of residues nearby the selected ligand. 
    '''
    make_sub_header('Getting maximum and minimum residues', out=self.logger)  
    maximum =  [key for key in self.res_dict if 
          all(self.res_dict[temp] <= self.res_dict[key]
          for temp in self.res_dict)]
            
    minimum =  [key for key in self.res_dict if 
            all(self.res_dict[temp] >= self.res_dict[key]
            for temp in self.res_dict)]

    
    print("Keys with minimum values are : " + str(minimum),"Keys with maximum values are : " + str(maximum))
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
  run_program(program_class=AnalyseGol)
# %%
