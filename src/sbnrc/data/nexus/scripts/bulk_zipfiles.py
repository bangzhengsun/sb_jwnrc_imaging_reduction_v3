import os, sys
import zipfile 


def zip_files(file_paths):
    """
    Given a list of file paths, create a .zip archive for each file.
    The archive will be named <original_filename>.zip and will contain the file.
    """
    for file_path in file_paths:
        if not os.path.isfile(file_path):
            print(f"Skipping: {file_path} (not a file)")
            continue

        # Build the zip filename by appending .zip
        zip_path = file_path + '.zip'

        # Create the zip and add the file
        with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            # arcname makes sure only the basename is stored inside the archive
            zf.write(file_path, arcname=os.path.basename(file_path))

        print(f"Created archive: {zip_path}")

if __name__ == "__main__": 
    files_to_zip = list(sys.argv[1:])
    zip_files(files_to_zip)