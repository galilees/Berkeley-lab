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

# TODO ******************
#results_dir = '/net/cci-filer3/home/galilees/pdb_survey_gol/''
results_dir ='/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/json_files/'
# TODO **********
# /net/cci-filer3/home/galilees/Berkeley-lab/run_galilee.py
script = '/net/anaconda/raid1/dorothee/14_frontiers_QR_restraints/galilee/run_galilee.py'
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
overwrite = False
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

  def run(self):
    '''
    '''
    #start_time = time.time()
    make_header('Running model %s' % self.pdb_code, out=self.logger)
    self.prepare_directory()
    self.initialize_json()
    #
    self.get_files_from_pdb_mirror()

    model = self.get_model_object(filename = self.json_data['pdb_file'])
    #model.overall_counts().show()

    #gol_dict_1 = {'nearby_res': [1,2,3,4], 'n_hbonds': 4}
    #self.json_data['GOL'] = {'sel_str1' : gol_dict_1}
    #self.save_json()
    #

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

#==============================================================================

if __name__ == '__main__':
  #
  from iotbx.cli_parser import run_program
  run_program(program_class=RunGenerate)
