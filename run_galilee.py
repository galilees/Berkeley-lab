from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import sys, os, json, datetime, time, six, shutil
import subprocess
import traceback
import iotbx.pdb
import mmtbx.model
from mmtbx import monomer_library
from elbow.quantum import electrons
from libtbx import easy_pickle
from libtbx import easy_run
from libtbx.utils import Sorry, null_out, to_str
from libtbx.str_utils import make_sub_header, make_header
from libtbx.command_line import easy_qsub
from iotbx import reflection_file_reader
from scitbx.array_family import flex


# TODO ******************
results_dir = '/net/cci-filer3/home/galilees/pdb_survey_gol/'
#results_dir ='/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/json_files/'
# TODO **********
script = '/net/cci-filer3/home/galilees/Berkeley-lab/run_galilee.py'
#script = '/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/run_galilee.py'
phenix_dir = "/net/cci-filer3/home/dcliebschner/phenix_svn/build/setpaths.csh"
#
pdb_folder = '/net/cci/pdb_mirror/pdb/'
mmcif_folder = '/net/cci/pdb_mirror/mmcif/'
#sfcif_folder = '/net/cci/pdb_mirror/structure_factors/'
# TODO
# file with list of pdb codes of models that contain GOL
#pickle_fn = '/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/scripts/pdb_codes_Xtallo_reso_range.pkl'
pickle_fn = ''

pdb_codes = ['1bg4']

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
mode = queue *one_cpu
  .type = choice
  .help = Run on the queue or one by one (one processor)
overwrite = True
  .type = bool
models
  .multiple=true
  {
  pdb_code = None
  .type = str
  }
'''

json_keys = [
'pdb_code',
'success',
'datum',
'GOL',
'HOH',
'pdb_file'
]

# =============================================================================

class RunGenerate(ProgramTemplate):

  datatypes = ['phil']
  master_phil_str = master_phil_str

#-----------------------------------------------------------------------------

  def validate(self):
    '''
    Simple validate
    '''
    make_sub_header('Validating inputs', out=self.logger)
    print('Running in mode: ', self.params.mode, file=self.logger)

    for model in self.params.models:
      code = model.pdb_code
      if len(code)!=4:
        raise Sorry("invalid PDB code for model")
      print(code, file=self.logger)

    # TODO
    #if not os.path.exists(pickle_fn):
    #  raise Sorry("pickle file does not exist: %s" % pickle_fn)
    #else:
    #  print('Using pickle file: ', pickle_fn, file=self.logger)

    print('Target directory: ' , results_dir, file=self.logger)

  #-----------------------------------------------------------------------------

  def run(self):
    '''
    Start queue jobs or run individually
    '''

    pdb_code_list = []
    if self.params.models:
      for model in self.params.models:
        pdb_code = model.pdb_code
        pdb_code_list.append(pdb_code)
    else:
      # TODO
      pdb_code_list = pdb_codes
      #if pdb_code != '1bg4': continue
      #pdb_code_list = easy_pickle.load(pickle_fn)
    #print(pdb_code_list[:5])
    #print(pdb_code_list)

    commands = list()
    n_jobs = 0
    #for pdb_code in pdb_code_list[:5]:
    for pdb_code in pdb_code_list:
      #
      if (self.params.mode == 'one_cpu'):
        obj = process_one_model(logger = self.logger,
                                pdb_code = pdb_code,
                                params = self.params)
        obj.run()
      if (self.params.mode == 'queue'):
        cmds = [
          'iotbx.python',
          script,
          'mode=one_cpu',
          'models=%s' % pdb_code
          ]
        cmd = " ".join(cmds)
        commands.append(cmd)
        #n_jobs += 1
        #if n_jobs==10000: break

    if (self.params.mode == 'queue'):
      pass
# TODO
# Uncomment this once the script is ready to be done on the queue
#
#      queue_log_dir = os.path.join(results_dir, 'queue_logs_' +
#        str(datetime.date.today()))
#      dir_name = queue_log_dir
#      if os.path.isdir(dir_name):
#        time_hour = str(datetime.datetime.now())[11:16]
#        time_hour = time_hour.replace(':', '_')
#        dir_name = queue_log_dir + '_' + time_hour
#      if (not os.path.isdir(dir_name)):
#        os.makedirs(dir_name)
#
#      for command in commands[:6]: print(command, file=self.logger)

#      easy_qsub.run(
#        phenix_source  = phenix_dir,
#        where          = queue_log_dir,
#        commands       = commands,
#        #qsub_cmd       = 'qsub -q all.q@%s -pe threaded 4' % machine,
#        qsub_cmd       = 'qsub -q all.q',
#        #qsub_cmd       = 'qsub -q all.q@morse',
#        js             = 20,
#        size_of_chunks = 1)

#==============================================================================

class process_one_model():
  datatypes = ['model']

  def __init__(self, logger, pdb_code, params):
    '''
      Initialize the class
    '''
    self.logger   = logger
    self.pdb_code = pdb_code
    self.params   = params

    self.success = True
    self.error_msg = ''
  
#-----------------------------------------------------------------------------  

  def initialize(self):
      '''
      inistialze data structures.
      '''
      #self.selection_dict = {}
      self.nearby_res_dict_gol = {}
      self.nearby_res_dict_hoh = {}
      self.selection_dict_gol = {}
      self.selection_dict_hoh = {}
      self.hbonds_dict_gol = {}
      self.hbonds_dict_hoh = {}
      self.nearby_residue_dict = {}
      self.res_dict = {}
      self.gol_hbonds_dict = {}
      self.hbonds_dict = {}
      self.json_data = {}
      # self.json_data['number of Hbonds involving GOL: '] = {}
      # self.json_data['nearby_res'] = {}

  #-----------------------------------------------------------------------------

  def run(self):
    '''
    '''
    self.success = True
    
    #
    self.initialize()
    #
    
    #start_time = time.time()
    make_header('Running model %s' % self.pdb_code, out=self.logger)
    self.prepare_directory()
    self.initialize_json()
    #
    self.get_files_from_pdb_mirror()

    self.model = self.get_model_object(filename = self.json_data['pdb_file'])
    #model.overall_counts().show()
    self.initialize()
    try:
      self.get_selection("GOL",self.selection_dict_gol)
    except Exception as e:
      self.success   = False
      print('failed to get GOL selection.\n' , file=self.logger)
      self.json_data["selection_strings"] = []
      self.save_json()

    try:
      self.get_selection("HOH",self.selection_dict_hoh)
    except Exception as e:
     self.success   = False
     print('failed to get HOH selection.\n' , file=self.logger)
     self.json_data["selection_strings"] = []
     self.save_json()
    
    try:  
      self.validate_gol(self.selection_dict_gol)
    except Exception as e:
      self.success   = False
      print('failed to validate gol selection.\n' , file=self.logger)
     

    try:  
      self.res_nearby_count(self.selection_dict_gol, self.nearby_res_dict_gol)
    except Exception as e:
      self.success   = False
      print('failed to get res nearby count for selection.\n' , file=self.logger)
      self.json_data["nearby_res"] = {}
      self.save_json()

    try:  
      self.res_nearby_count(self.selection_dict_hoh, self.nearby_res_dict_hoh)
    except Exception as e:
      self.success   = False
      print('failed to get res nearby count for selection.\n' , file=self.logger)
      self.json_data["nearby_res"] = {}
      self.save_json()

    
    try:
      self.get_hbonds(self.selection_dict_gol)
    except Exception as e:
      self.success   = False
      print('failed to get hbonds for selection.\n' , file=self.logger)
      self.json_data['number of Hbonds involving GOL: '] = {}
      self.save_json()

    try:
      self.get_hbonds(self.selection_dict_hoh)
    except Exception as e:
      self.success   = False
      print('failed to get hbonds for selection.\n' , file=self.logger)
      self.json_data['number of Hbonds involving GOL: '] = {}
      self.save_json()
       

    # pdb_hierarchy = self.model.get_hierarchy()
    # sec_str_from_pdb_file = self.model.get_ss_annotation()

    # # get secodary structure annotation vector from HELIX/SHEET records (file header)
    # print('Running secondary structure annotation...')
    # v1 = self.get_ss(hierarchy             = pdb_hierarchy,
    # sec_str_from_pdb_file = sec_str_from_pdb_file)
    
    
  
    # gol_dict_1 = {json.dumps(self.nearby_residue_dict.items()):json.dumps(self.hbonds_dict.items())

    gol_dict_json = {json.dumps({"nearby res" : self.nearby_res_dict_gol.values()}) : json.dumps({'number of Hbonds involving GOL: ': self.hbonds_dict_gol.values()})}
 
   
    self.json_data['GOL'] = {json.dumps(self.selection_dict_gol.keys()):gol_dict_json}

    hoh_dict_json = {json.dumps({"nearby res" : self.nearby_res_dict_hoh.values()}) : json.dumps({'number of Hbonds involving HOH: ' : self.hbonds_dict_hoh.values()})}
 
   
    self.json_data['HOH'] = {json.dumps(self.selection_dict_hoh.keys()):hoh_dict_json}
    self.save_json()

    # key = frozenset(self.nearby_residue_dict.items())
    # value = frozenset(self.hbonds_dict.items())
    # gol_dict_1 = (key,value)
    # key2 = self.selection_dict.keys()
    # self.json_data['GOL'] = {key2 : gol_dict_1}
    # self.save_json()
    
    
    # gol_dict_1 = {'nearby_res': [1,2,3,4], 'n_hbonds': 4}
    #self.json_data['GOL'] = {'sel_str1' : gol_dict_1}
    # self.save_json()

    

  #-----------------------------------------------------------------------------

  def get_model_object(self, filename):
    '''
      Get model object without restraints.
    '''
    model = None
    make_header('Getting model object', out=self.logger)
    try:
      pdb_inp = iotbx.pdb.input(file_name = filename)
    except Exception as e:
      msg = traceback.format_exc()
      print(msg, file=self.logger)
      self.success = False
      self.write_log()
      return
    try:
      model = mmtbx.model.manager(model_input = pdb_inp)
    except Exception as e:
      msg = traceback.format_exc()
      print(msg, file=self.logger)
      self.success = False
      self.write_log()
      return
    print('...finished', file=self.logger)
    return model

  #-----------------------------------------------------------------------------

  def get_files_from_pdb_mirror(self):
    '''
    Fetch model and data from the pdb mirror
    '''
    make_sub_header('Get model and data from mirror', out=self.logger)
    #model_fn = self.model_fn_mirror
    model_fn = self.get_mirror_fn(filedir  = pdb_folder)
    print(model_fn)
    pdb_fn = os.path.join(self.dest_dir, self.pdb_code +'.pdb')
    #mtz_fn = os.path.join(self.dest_dir, self.pdb_code +'.mtz')
    #basename, extension = os.path.splitext(os.path.basename(model_fn))
    if not os.path.exists(pdb_fn):
      model_fn_dest_dir = os.path.join(self.dest_dir, self.pdb_code + '.pdb.gz')
      print(model_fn_dest_dir)
      print('Copying and unpacking ', model_fn, file=self.logger)
      shutil.copyfile(model_fn, model_fn_dest_dir)
      cmds = ["gunzip", model_fn_dest_dir]
      print(" ".join(cmds), file=self.logger)
      try:
        r = subprocess.check_output(cmds)
      except subprocess.CalledProcessError as e:
        msg = e.output
    else:
      print('Found file %s' % pdb_fn, file=self.logger)
    #
#    if not os.path.exists(mtz_fn):
#      sfcif_fn = self.get_mirror_fn(filedir = sfcif_folder)
#      sfcif_fn_dest_dir = os.path.join(self.dest_dir, self.pdb_code + '-sf.cif')
#      sfcif_fn_dest_dir_gz = sfcif_fn_dest_dir + '.gz'
#      print('Copying and unpacking ', sfcif_fn, file=self.logger)
#      shutil.copyfile(sfcif_fn, sfcif_fn_dest_dir_gz)
#      cmds = ["gunzip", sfcif_fn_dest_dir_gz]
#      print(" ".join(cmds), file=self.logger)
#      try:
#        r = subprocess.check_output(cmds)
#      except subprocess.CalledProcessError as e:
#        msg = e.output
#      cmds = ["phenix.cif_as_mtz",
#        sfcif_fn_dest_dir,
#        "--merge",
#        "--output_file_name=%s" % self.pdb_code+'.mtz']
#      print(" ".join(cmds), file=self.logger)
#      try:
#        r = subprocess.check_output(cmds)
#      except subprocess.CalledProcessError as e:
#        msg = e.output
#    else:
#      print('Found file %s' % mtz_fn, file=self.logger)
    #
    if os.path.isfile(pdb_fn): self.json_data['pdb_file'] = pdb_fn
    else: self.success = False
#    if os.path.isfile(mtz_fn): self.json_data['mtz_file'] = mtz_fn
#    else: self.success = False
    self.save_json()

  #-----------------------------------------------------------------------------

  def get_mirror_fn(self, filedir):
    '''
    '''
    model_file = None
    fo=open(filedir+"/INDEX","r")
    for pdb_file_ in fo.readlines():
      pdb_file_ = pdb_file_.strip()
      if(pdb_file_.count(self.pdb_code)):
        model_file = filedir+pdb_file_
        break
    fo.close()
    return model_file

  #----------------------------------------------------------------------------

  def save_json(self):
    '''
    Save results in json file
    '''
    print('Saving json file.', file=self.logger)
    self.json_data['success'] = self.success
    self.json_data['error_msg'] = self.error_msg
    with open(self.json_fn, 'w') as fp:
      json.dump(self.json_data, fp, sort_keys=True, indent=4)

  #-----------------------------------------------------------------------------

  def initialize_json(self):
    '''
    Initialize json file
    '''
    self.json_fn = os.path.join(self.dest_dir, self.pdb_code + '.json')
    if (os.path.isfile(self.json_fn)):
      print('Opened file', os.path.basename(self.json_fn), file=self.logger)
      with open(self.json_fn, 'r') as fp:
        self.json_data = json.load(fp)
        for key in json_keys:
          if key not in self.json_data:
            self.json_data[key] = None
    else:
      self.json_data = dict()
      for key in json_keys:
        self.json_data[key] = None
    #
    self.json_data['success'] = True
    self.json_data['pdb_code'] = self.pdb_code
    self.json_data['datum']    = str(datetime.date.today())
    self.save_json()
    print('Initialized json', file=self.logger)

  #----------------------------------------------------------------------------

  def prepare_directory(self):
    '''
    Prepare directory for each job
#    '''
    self.dest_dir = results_dir
#    make_sub_header('Preparing the result directory', out=self.logger)
#    dest_dir = os.path.join(results_dir, self.pdb_code)
#    # Make new directory if it does not exist
#    if (not os.path.isdir(dest_dir)):
#      os.makedirs(dest_dir)
#      print('Making new directory', dest_dir, file=self.logger)
#    # If overwrite=True and directory exists, delete all content in the directory
#    elif (self.params.overwrite):
#      print('Deleting folder content', file=self.logger)
#      for one_file in os.listdir(dest_dir):
#          file_path = os.path.join(dest_dir, one_file)
#          try:
#              if os.path.isfile(file_path):
#                  os.unlink(file_path)
#          except Exception as e:
#              print('Deleting folder content unsuccessful.', file=self.logger)
#              print(e, file=self.logger)
#    else:
#      print('Using existing directory: ', dest_dir, file=self.logger)
#    os.chdir(dest_dir)
#    self.dest_dir = dest_dir

 #----------------------------------------------------------------------------
  def get_ss(self,
           hierarchy,
           sec_str_from_pdb_file=None,
           method="ksdssp",
           use_recs=False):
    if(use_recs): params = None
    else:
      params = mmtbx.secondary_structure.manager.get_default_ss_params()
      params.secondary_structure.protein.search_method=method
      params = params.secondary_structure
    ssm = mmtbx.secondary_structure.manager(
      pdb_hierarchy         = hierarchy,
      sec_str_from_pdb_file = sec_str_from_pdb_file,
      params                = params,
      log                   = null_out())
    #print(ssm.records_for_pdb_file())
    alpha = ssm.helix_selection()
    beta  = ssm.beta_selection()
    #print(list(alpha))

    #print(dir(ssm))
    #STOP()

    assert alpha.size() == beta.size() == hierarchy.atoms().size()
    annotation_vector = flex.double(hierarchy.atoms().size(), 0)
    annotation_vector.set_selected(alpha, 1)
    annotation_vector.set_selected(beta, 2)

    #records = ssm.records_for_pdb_file()
    #print(records)
    sel_str = 'chain A and resname GOL and resseq  630'
    near_res = "residues_within(5,%s)"%sel_str
    selection_bool2 = self.model.selection(near_res)
    m2 = self.model.select(selection_bool2)

    m2.set_ss_annotation(ann = sec_str_from_pdb_file)
    print(m2.model_as_pdb())
    # sec_str_from_m2 = m2.get_ss_annotation()
    # print(sec_str_from_m2)
    params = mmtbx.secondary_structure.manager.get_default_ss_params()
    params.secondary_structure.protein.search_method="ksdssp"
    params = params.secondary_structure
    params = None
    ssm = mmtbx.secondary_structure.manager(
      pdb_hierarchy         = m2.get_hierarchy(),
      sec_str_from_pdb_file = sec_str_from_pdb_file,
      params                = params,
      log                   = null_out())
   
    print(ssm.helix_selection().count(True))
    print(len(ssm.helix_selection()))
    print(ssm.get_helix_types)
    print(ssm.find_approximate_helices())

    return annotation_vector



 #-----------------------------------------------------------------------------

  def print_hbond_table(self, model, gol_hbonds_dict):
    '''
    takes in the model and a dictionary of hbonds and prints the table associated with the hbonds
    '''
    result_str = '{:<18} : {:5d}'
    # print table with all H-bonds
    title1 = ['donor', 'acceptor', 'distance', 'angle']
    title1_str = '{:^33}|{:^16}|{:^21}|{:^14}|'
    print('\n' + title1_str.format(*title1))
    title2 =  ['X', 'H', 'A','H...A','X...A',
              'X-H...A', 'symop']
    title2_str = '{:^16}|{:^16}|{:^16}|{:^10}|{:^10}|{:^14}|{:^15}|'
    print(title2_str.format(*title2))
    table_str = '{:>16}|{:>16}|{:^16}|{:^10.2f}|{:^10.2f}|{:^14.2f}|{:^15}|'
    print('-'*99)
    ###########
    atoms = model.get_atoms()
    for iseq_tuple, record in gol_hbonds_dict.items():
      iseq_x, iseq_h, iseq_a = iseq_tuple
      if record[4] is not None:
        symop = record[4]
      else: symop = ''
      x_id_str = atoms[iseq_x].id_str().replace('pdb=','').replace('"','')
      h_id_str = atoms[iseq_h].id_str().replace('pdb=','').replace('"','')
      a_id_str = atoms[iseq_a].id_str().replace('pdb=','').replace('"','')
      line = [x_id_str, h_id_str, a_id_str, round(record[0], 2),
        round(record[1], 2), round(record[2], 2), symop]
      print(table_str.format(*line))
    print('-'*99)

  #-----------------------------------------------------------------------------

  def get_hbonds(self,hbonds_dict):
    '''
    Searches within a radius of 5 angstroms for hydrogen bonds to glycerol and returns a list of tuples containing the isequences of hbond partners 
    '''
  
    make_sub_header('H-bonds', out=self.logger)
   
    hbonds_list=[]
    iselection_dict = {}
    _i = 0
    for sel_str in self.selection_dict.keys():
      _i += 1
      #if _i > 1: break
      print('Now looking at ', sel_str)
      near_res = 'residues_within(5,%s)'%sel_str
      selection_bool2 = self.model.selection(near_res)
      m2 = self.model.select(selection_bool2)
      m2.set_log(log = null_out())
      m2.process(make_restraints=True)

      pnps = pnp.manager(model = m2)
    
      hbonds = pnps.get_hbonds()
      #hbonds.show(log=sys.stdout) # hide
      print("tuple keys to hbonds table: ", hbonds._hbonds_dict.keys()) #Show tuple index to table
    
      # better name? gol_iseq_numbers
      iselection_list = list(m2.iselection(sel_str))

      gol_hbonds_dict = {}
      for iseq_tuple in list(hbonds._hbonds_dict.keys()):
        for iseq in iseq_tuple:
          if iseq in iselection_list:
            gol_iseq_tuple = iseq_tuple
            hbonds_dict_value = hbonds._hbonds_dict.get(gol_iseq_tuple)
            gol_hbonds_dict[gol_iseq_tuple] = hbonds_dict_value
            if gol_iseq_tuple not in hbonds_list:
              hbonds_list.append(gol_iseq_tuple)
          else:
            break

      self.print_hbond_table(model           = m2,
                             gol_hbonds_dict = gol_hbonds_dict)
      hbonds_dict[sel_str] = len(gol_hbonds_dict.keys())

      # self.json_data['number of Hbonds involving GOL: '] = self.hbonds_dict
      # self.save_json()
      print('number of Hbonds involving GOL: ', len(gol_hbonds_dict.keys()))
     
    
 #-----------------------------------------------------------------------------  

  def get_selection(self, residue, selection_dict):
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
             
              selection_dict[sel_str] = iselection

              
      # self.json_data['selection_strings'] = list(self.selection_dict.keys())
      # self.save_json()
    
    for sel_str in self.selection_dict.keys():
      print(sel_str)
    
              
#----------------------------------------------------------------------------
  def validate_gol(selection_dict,self, b_max = 100,occ_min = .2):
    '''
    b_max: int
    occ_min: int
    removes strutures from the selection dictionary that do no meet criteria. The b factor is out of range if it exceeds b_max and the occupancy is out of range if it is less than occ_min. 
    Also checks that all occupancies within the molecule are the same.
    '''

    make_sub_header('curate gol selection', out=self.logger)
    
    bad_selection=[]
   
    for sel_str in selection_dict.keys():
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

  def res_nearby_count(self,selection_dict,nearby_residue_dict):
    '''
   counts the residues nearby each selection. Returns a dictionary with the name of the residue as the key and its count as the value. 
    '''
    make_sub_header('count of residues nearby each selection', out=self.logger)
    
    
    for sel_str in selection_dict.keys():
     
      near_res = "residues_within(5,%s)"%sel_str
      selection_bool2 = self.model.selection(near_res)
      m = self.model.select(selection_bool2) 
      ph = m.get_hierarchy()
      res_nearby = ph.overall_counts().resnames
      #print(res_nearby, file=self.logger)
      nearby_residue_dict[sel_str] = res_nearby

    # self.json_data['nearby_res'] = self.nearby_residue_dict
    # self.save_json()
    # return res_nearby


if __name__ == '__main__':
  #
  from iotbx.cli_parser import run_program
  run_program(program_class=RunGenerate)
 
 
 
 
 
  

#==============================================================================