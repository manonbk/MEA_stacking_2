from Bio import AlignIO

# Define the path to the test files
import sys
import os

# Remonter d'un niveau pour atteindre le dossier parent (MEA_STACKING)
parent_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
sys.path.append(parent_dir)

# Définir le chemin du dossier 'data'
folder_path = os.path.join(parent_dir, "data")

print("Folder path:", folder_path)



def read_alignment(file_path, file_format):
    """
    Reads an RNA alignment from a given file in Clustal or Stockholm format.

    :param file_path: Path to the alignment file.
    :param file_format: Format of the alignment ('clustal' or 'stockholm').
    :return: Biopython MultipleSeqAlignment object.
    """
    try:
        alignment = AlignIO.read(file_path, file_format)
        return alignment
    except Exception as e:
        print(f"Error reading alignment file: {e}")
        return None



def detect_format(filename):
    """
    Detectes file format according to the extension
    """
    if filename.endswith(".stockholm.txt") or filename.endswith(".stk"):
        return "stockholm"
    elif filename.endswith(".aln") or "clustalw" in filename:
        return "clustal"
    elif filename.endswith(".fa"):
        return "fasta"
    else:
        return None  # Format inconnu


def import_test_files():
    """
    Imports all the test alignements in the test_data folder.
    Returns a list of alignement, each element in the list corresponds to a file, 
    and is the list of the aligned sequences from the file.
    """
    test_aligned_sequences={}
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        
        # Détecter le format
        file_format = detect_format(filename)
        
        if file_format:  # Vérifier que le format est supporté
            try:
                alignment = AlignIO.read(file_path, file_format)
                sequences = [str(record.seq) for record in alignment]
                test_aligned_sequences[filename]=sequences
                print(f"✅ File {filename} loaded with success ({file_format})")
            except Exception as e:
                print(f"❌ Error while reading {filename} : {e}")
        else:
            print(f"⚠ Error of file format {filename}, ignored.")

    return test_aligned_sequences


#print(import_test_files())
