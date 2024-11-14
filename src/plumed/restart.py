# TODO: this is a mess, needs to be cleaned up


# print("Getting CVs from previous plumed.dat")
# # assume it's in the same folder as restart_rfile
# restart_rfile_path = os.path.dirname(restart_rfile)
# plumed_file = os.path.join(restart_rfile_path, f'{filename}_plumed.dat')
# assert os.path.exists(plumed_file), f"File {plumed_file} does not exist"
# shutil.copy(plumed_file, f"{output_dir}/{filename}_plumed.dat")

# # assume fixed.pdb also in the same folder
# fixed_pdb = os.path.join(restart_rfile_path, f'{filename}_fixed.pdb')
# assert os.path.exists(fixed_pdb), f"File {fixed_pdb} does not exist"
# # if it is, copy it to the tmp folder we deal with
# shutil.copy(fixed_pdb, f'{output_dir}/{filename}_fixed.pdb')

# # do the same copying for _solvated.pdb
# solvated_pdb = os.path.join(restart_rfile_path, f'{filename}_solvated.pdb')
# assert os.path.exists(solvated_pdb), f"File {solvated_pdb} does not exist"
# shutil.copy(solvated_pdb, f'{output_dir}/{filename}_solvated.pdb')

# # TODO: do I need to copy the kernels? the colvar? anything like that?

# restart_checkpoint = os.path.join(restart_rfile_path, f'{filename}.chk')
# assert os.path.exists(restart_checkpoint), f"File {restart_checkpoint} does not exist"

# def extract_contact_pairs_str(plumed_file):
#     import re
#     with open(plumed_file, 'r') as f:
#         content = f.read()
#     # Find all matches
#     matches = re.findall(r'ATOMS\d+=\d+,\d+', content)

#     # Join matches into a single string with newlines
#     result_string = '\n\t'.join(matches)
#     result_string = f"	{result_string}"  # Add leading tab for formatting
    
#     return result_string

# contact_pairs_str = extract_contact_pairs_str(plumed_file)
