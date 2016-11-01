import os
import sys
import subprocess
import lt
for f in lt.files_in_dir(sys.argv[-1]):
    subprocess.call('~/hbplus/hbplus/hbplus',f,"-f 'SEP.txt'")
