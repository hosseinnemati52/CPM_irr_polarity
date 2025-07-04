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


def create_zip_from_folder(folder_path, zip_file_name):
    try:
        with zipfile.ZipFile(zip_file_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(folder_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    archive_name = os.path.relpath(file_path, folder_path)
                    zipf.write(file_path, archive_name)

        print(f"Successfully created {zip_file_name} from {folder_path}")
    except Exception as e:
        print(f"Error: {e}")

def clear_folder(folder_path):
    try:
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
        print(f"Cleared the contents of {folder_path}")
    except Exception as e:
        print(f"Error: {e}")

# zipFileName = 'data.zip'
# folderName = 'data'

zipFileName = sys.argv[1]
folderName = sys.argv[2]

# with zipfile.ZipFile(zipFileName, 'w', zipfile.ZIP_DEFLATED) as zipf:
#     zipdir(folderName+'/', zipf)

create_zip_from_folder(folderName, zipFileName)

with zipfile.ZipFile(zipFileName, 'r') as zip_ref:
    
    # check the integrity of the zipfile
    result = zip_ref.testzip()
    
    # print the result of the integrity check
    if result is not None:
        print(f"Error: zip file "+zipFileName+" is corrupted")
        sys.exit()
    else:
        # print("Zip file is good, extracting...")
        # zip_ref.extractall('')
        # print("Zip file extracted!")
        print("Zip file "+zipFileName+" is good! Folder content is deleted!")
        # shutil.rmtree('data', ignore_errors=False, onerror=None)
        clear_folder('data')
