import os
import subprocess


class GeneratePeptides:
    """
    Attributes
    ==========
    sequences : list
        A list of sequences to create the cyclic structure with.

    Methods
    =======
    Generate(self, sequences)
        Use rosetta to generate the structure of a cylcic peptide.

    Hidden Methods
    ==============

    """


    def __init__(self, sequences):
        """
        Check input sequences, transform to suitable format and create input files.

        Parameters
        ==========
        sequences : str/list
            Sequences to create the cyclic peptides.
        """

        self.sequences = sequences

    def _letter_to_word_aa(self, sequence):
        """
        Change fromat from a string of amino acid letters to the correct format.

        Parameters
        ==========
        sequence : str
            string of letters of amino acids.
        """

        aa_dict = {
            'A':'ALA',
            'C':'CYS',
            'D':'ASP',
            'E':'GLU',
            'F':'PHE',
            'G':'GLY',
            'H':'HIS',
            'I':'ILE',
            'K':'LYS',
            'L':'LEU',
            'M':'MET',
            'N':'ASN',
            'P':'PRO',
            'Q':'GLN',
            'R':'ARG',
            'S':'SER',
            'T':'THR',
            'V':'VAL',
            'W':'TRP',
            'Y':'TYR'
        }
        aa_list = []
        for aa in sequence:
            aa_list.append(aa_dict[aa])
        sequence_formatted = ' '.join(aa_list)
        sequence_formatted += '\n'

        return sequence_formatted

    def _create_file(self, sequence, out_dir, name):
        """
        Create the input files with the correct sequence_format.

        Parameters
        ==========
        sequence : str
            Sequence to which a file will be created
        out_dir : str
            Directory where output_files are placed
        name : str
            Name of the output_files
        """

        if out_dir not in os.listdir():
            os.system('mkdir {}'.format(out_dir))

        input_sequence = self._letter_to_word_aa(sequence)
        with open(out_dir+'/'+name, 'w') as f:
            f.write(input_sequence)

    def Generate(self, output_dir, iterations=100000):
        """
        Run the rosetta software to create cyclic peptides structures.

        Parameters
        ==========
        input_file : str
            file from where sequence is taken
        output_dir : str
            name of the directory where the structures will be written
        iterations : int
            number of attempts to create the structure
        """

        def _simplify_pdb(input_file, output_file, conect=True):
            """
            Remove not important lines from the created PDB and add connect line to close the cycle.
            """
            atom_lines = []
            connect_lines = []
            other_lines = []

            with open(input_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        atom_lines.append(line)
                    elif line.startswith('CONECT'):
                        connect_lines.append(line)
                    elif line.startswith('LINK'):
                        other_lines.append(line)

            if conect:
                residues_to_conect = (other_lines[0].split()[4], other_lines[0].split()[8])
                with open(input_file, 'r') as f:
                    for line in f:
                        if line.startswith('ATOM'):
                            #print(line.split())
                            if line.split()[2] == 'N' and line.split()[5] == str(residues_to_conect[0]):
                                atom_1 = line.split()[1]
                            elif line.split()[2] == 'C' and line.split()[5] == str(residues_to_conect[1]):
                                atom_2 = line.split()[1]

            with open(output_file, 'w') as f:
                f.write(''.join(atom_lines))
                f.write('TER\n')  # Add newline before connect lines
                f.write(''.join(connect_lines))
                conect_spaces = ''.join((5-len(atom_1))*[' '])
                f.write('CONECT'+conect_spaces+atom_1+'  '+atom_2+'\n')
                f.write('END\n')

        if output_dir not in os.listdir():
            os.system('mkdir {}'.format(output_dir))
        count = 0
        for s in self.sequences:
            attempts = iterations
            count += 1
            name = 'test'+str(count)
            self._create_file(s, 'input_files', name+'.txt')

            input_file = 'input_files/'+name+'.txt'
            output_file = output_dir+'/raw_'+name

            print('Calculating structure for: ')
            os.system('cat {}'.format(input_file))

            command = 'simple_cycpep_predict.default.linuxgccrelease'
            command += ' -cyclic_peptide:sequence_file '+input_file
            command += ' -cyclic_peptide:genkic_closure_attempts '+str(attempts)
            command += ' -mute '+ 'all'
            command += ' -unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication'
            command += ' -out:file:o '+output_file

            try:
                output = subprocess.check_output(
                            command, stderr=subprocess.STDOUT, shell=True,
                            universal_newlines=True)
                print(output)

            except subprocess.CalledProcessError as exc:
                print(exc.output, exc.returncode)
                raise Exception("Problem with input")

            if os.path.isfile(output_file+'.pdb') == False:
                print('No structure generated for this peptide')

            elif os.path.isfile(output_file) == True:
                print('done')

        ## Add last connect line to close the cycles
        files = [x for x in os.listdir(output_dir) if x.startswith('raw')]


        for f in files:
            print(f)
            _simplify_pdb(output_dir+'/'+f, output_dir+'/'+f[4:])
            os.system('rm {}'.format(output_dir+'/'+f))
