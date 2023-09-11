from classy_szfast.classy_szfast import Class_szfast
# from classy_szfast.classy_sz import classy_sz
from .config import *
from .utils import *
from .cosmopower import *
from .pks_and_sigmas import *

from pathlib import Path
path = str(Path(__file__).parent.absolute())
import sys
path = path+'/../../../class_sz_auxiliary_files/custom_profiles'
sys.path.insert(0,path)
# print(path)
from custom_profiles import *
