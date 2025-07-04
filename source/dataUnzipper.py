#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 10:37:14 2023

@author: hossein
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

# with zipfile.ZipFile('data.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
#     zipdir('data/', zipf)


with zipfile.ZipFile('data.zip', 'r') as zip_ref:
    
    # check the integrity of the zipfile
    result = zip_ref.testzip()
    
    # print the result of the integrity check
    if result is not None:
        print(f"Error: zip file is corrupted")
        sys.exit()
    else:
        print("Zip file is good, extracting...")
        zip_ref.extractall('')
        print("Zip file extracted!")
        # print("Zip file is good! Folder is deleted!")        
        # shutil.rmtree('data', ignore_errors=False, onerror=None)

switch = sys.argv[1]

if switch=='all':

    data_index = 1
    while 1:
        zipFileName = 'data_'+str(data_index)+'.zip'
        try:

            with zipfile.ZipFile(zipFileName, 'r') as zip_ref:
        
                # check the integrity of the zipfile
                result = zip_ref.testzip()
                
                # print the result of the integrity check
                if result is not None:
                    print(f"Error: zip file "+zipFileName+" is corrupted")
                    sys.exit()
                else:
                    print("Zip file "+zipFileName+" is good, extracting...")
                    zip_ref.extractall('data')
                    print("Zip file "+zipFileName+" extracted!")

                data_index += 1

        except:
            break

elif switch =='data':
    pass

else:
    print("Wrong argvalue for dataUnzipper")
    exit()
