import pandas as pd


class AlignmentInformation(object):
    """This class uses a class Reference instance to initialize a count_dictionary that stores alignment nucleotide count
    information. count_dictionary only accounts information from exon positions.
    By reading a class Read instance, it updates nucleotide count information to the count_dictionary.
    Finally, method create_proportion_dictionary creates a proportion_dictionary,
    that calculate nucleotide proportion per position, using count information.
    ...

    Attributes
    ----------
    reference : instance from class Reference
        A class Reference instance
    count_dictionary : dict
        A dictionary of dictionaries accessed by position in alignment keys. Each position key contains five keys,
        one for each possible nucleotide in the alignment: A, T, G, C and . (gap)
        Every nucleotide key relates with counter (float) values that are initialized as 0.0
    proportion_dictionary : dict
        A dictionary of dictionaries accessed by position in alignment keys. Each position key contains five keys,
        one for each possible nucleotide in the alignment. Every nucleotide key has a proportion (float) value.
    reference_alleles : list
        List of names of alleles in the reference alignment
    processed_reads_dictionary : dict
        Dictionary that stores already processed reads as keys. Values are length based positions ranges of the reads
    count_dataframe : pandas dataframe
        count_dictionary converted into a pandas dataframe, for manual inspection
    proportion_dataframe : pandas dataframe
        proportion_dictionary converted into a pandas dataframe, for manual inspection

    Methods
    -------
    update_count_dictionary
        Takes a class Read instance as an argument and updates the count_dictionary
        It considers if the given read is overlapping with its paired-end, to avoid repetition
    create_proportion_dictionary
        Creates a proportion_dictionary out of count_dictionary. Per position in count_dictionary, it calculates the
        proportions of every nucleotide.
    print_dictionaries_as_dataframe
        Prints count_dictionary and/or proportion_dictionary as pandas dataframes for manual inspection

    """

    def __init__(self, a_reference_instance):
        """Count_dictionary is initialized using exon indexes from a class Reference instance. Only exon positions
        are included

        Parameters
        ----------
        a_reference_instance
            Instance inherited from class Reference
        """

        self.reference = a_reference_instance
        self.count_dictionary = {}
        self.proportion_dictionary = {}
        self.count_dataframe = None
        self.proportion_dataframe = None
        for exon_range in self.reference.exons_index_list:
            for i in range(exon_range[0], exon_range[1]):
                self.count_dictionary[i] = {"A": 0.0, "T": 0.0, "G": 0.0, "C": 0.0, ".": 0.0}

        self.reference_alleles = []
        with open(self.reference.file_name) as alignment_reference:
            for allele in alignment_reference.readlines()[12:-1]:  # First sequence in reference starts in line 12
                allele = allele.replace(" ", "")
                allele_id = allele.split("\t")[0]
                self.reference_alleles.append(allele_id)
        alignment_reference.close()
        self.processed_reads_dictionary = {}

    def update_count_dictionary(self, a_read_instance):
        """Takes a class Read instance, check whether this read's pair was already processed
        Updates count_dictionary, skipping already accounted information from paired-end reads

        Parameters
        ----------
        a_read_instance
            Instance inherited from class Read

        Raises
        ------
        If a nucleotide is defined as "N", no information is updated to the count_dictionary in that position
        """

        position = a_read_instance.leftmost_position
        if a_read_instance.id in self.processed_reads_dictionary:
            paired_read_range = set(range(a_read_instance.leftmost_position, a_read_instance.rightmost_position))
            paired_reads_overlapping_positions = paired_read_range.intersection(
                range(self.processed_reads_dictionary[a_read_instance.id][0],
                      self.processed_reads_dictionary[a_read_instance.id][1]))
        else:
            self.processed_reads_dictionary[a_read_instance.id] = [a_read_instance.leftmost_position,
                                                                   a_read_instance.rightmost_position]
            paired_reads_overlapping_positions = []
        for nucleotide in a_read_instance.aligned_sequence:
            if nucleotide == "N":
                pass
            elif position in self.count_dictionary and position not in paired_reads_overlapping_positions:
                self.count_dictionary[position][nucleotide] += 1
            position += 1

    def create_proportion_dictionary(self):
        """Creates a proportion_dictionary out of count_dictionary. Per position in dictionary, it calculates the
        proportions of every nucleotide. Proportions are rounded up to two decimals.

        """

        for position in self.count_dictionary:
            nucleotide_count = sum(self.count_dictionary[position].values())
            self.proportion_dictionary[position] = {}
            for nucleotide in self.count_dictionary[position]:
                self.proportion_dictionary[position][nucleotide] = round(self.count_dictionary[position][
                                                                           nucleotide] / nucleotide_count * 100, 2)

    def print_dictionaries_as_dataframe(self, print_count_dictionary=True, print_proportion_dictionary=False):
        """Creates Pandas DataFrame from count_dictionary and/or proportion_dictionary
        Prints dataframe or dataframes for manual inspection

        Parameters
        ----------
        print_count_dictionary : boolean
            Indicates whether count_dictionary is to be printed as dataframe
        print_proportion_dictionary : boolean
            Indicates whether proportion_dictionary is to be printed as dataframe

        Raises
        ------
        If a proportion_dictionary was not created using the method create_proportion_dictionary, a warning is shown
        """

        if print_count_dictionary is True:
            self.count_dataframe = pd.DataFrame(self.count_dictionary)
            pd.set_option('display.max_columns', None)
            print(self.count_dataframe)
        if print_proportion_dictionary is True:
            if self.proportion_dictionary == {}:
                return "WARNING: a proportion dictionary is not created"
            else:
                self.proportion_dataframe = pd.DataFrame(self.proportion_dictionary)
                pd.set_option('display.max_columns', None)
                print(self.proportion_dataframe)
