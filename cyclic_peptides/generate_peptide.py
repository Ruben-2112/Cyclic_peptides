import os


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

    def Generate(self, sequences):
        """
        Run the rosetta software to create cyclic peptides structures.

        Parameters
        ==========

        """

        print(sequences)
