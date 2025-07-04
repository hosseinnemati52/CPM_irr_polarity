# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 10:37:39 2023

@author: Nemat002
"""

import os
import zipfile
import sys
import shutil
    
def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file), 
                       os.path.relpath(os.path.join(root, file), 
                                       os.path.join(path, '..')))

def zipdir_general(path, ziph):
    # ziph is zipfile handle
    flag = 0
    for root, dirs, files in os.walk(path):
        for file in files:
            flag = 1
            ziph.write(os.path.join(root, file), 
                       os.path.relpath(os.path.join(root, file), 
                                       os.path.join(path, '..')))
    if flag==0:
        ziph.writestr("data/dummy.txt", "This is a dummy file. :)")

with zipfile.ZipFile('data.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
    # zipdir('data/', zipf)
    zipdir_general('data/', zipf)
    # pass


with zipfile.ZipFile('data.zip', 'r') as zip_ref:
    
    # check the integrity of the zipfile
    result = zip_ref.testzip()
    
    # print the result of the integrity check
    if result is not None:
        print(f"Error: zip file is corrupted")
        sys.exit()
    else:
        # print("Zip file is good, extracting...")
        # zip_ref.extractall('')
        # print("Zip file extracted!")
        print("Zip file is good! Folder is deleted!")        
        shutil.rmtree('data', ignore_errors=False, onerror=None)
