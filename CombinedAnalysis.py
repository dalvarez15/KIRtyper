import itertools as it
import numpy as np


class CombinedAnalysis(object):
    """This class perform Combined Analysis using a class ProgressiveAnalysis instance result_alleles_dictionary attribute,
    a class AlignmentInfo instance and reference file.

    ...

    Attributes
    ----------
    all_genotype_combinations : list of lists of tuples
        All possible genotype combinations (every list inside the final list)
        of 2 alleles (every tuple inside every list) per detected KIR gene
    coding_sequences_per_position : dictionary of dictionaries
        Per exon position, it stores nucleotide information from every result allele from class Progressive Analysis
    genotype_combinations_per_position : list of subsets of all_genotype_combinations
        Every item in the list is a subset of all_genotype_combinations that matched every analyzed alignment position
    typing_result : list of lists of tuples
        Subset of all_genotype_combinations that were found to match every position in the class AlignmentInfo instance
    result_alleles_dictionary : dict
        Dictionary of lists, one per detected KIR gene, storing names of result alleles from Progressive Analysis

    Methods
    -------
    get_result_alleles_coding_sequences
        Gets coding sequences of alleles in a class ProgressiveAnalysis result_alleles_dictionary attribute
    find_matching_combinations_per_position
        Calls methods is_discriminant() and check_genotype_combinations()
        to find perfect matching combinations per position
    is_discriminant
        It determines whether a position in the alignment is discriminant
    check_genotype_combinations
        It determines which genotype combinations match the alignment information, per analyzed position
    get_combined_results
        Find all genotype combinations that match all discriminant possible and get final_result
    """

    def __init__(self, progressive_analysis_instance, alignment_information_instance):
        """Calls combined analysis methods

        Parameters
        ----------
        progressive_analysis_instance : class ProgressiveAnalysis instance
            Instance inherited from class ProgressiveAnalysis
        alignment_information_instance : class AlignmentInfo instance
            Instance inherited from class ProgressiveAnalysis
        """
        self.all_genotype_combinations = []
        self.coding_sequences_per_position = {}
        self.genotype_combinations_per_position = []
        self.typing_result = None
        self.result_alleles_dictionary = {}

        self.get_result_alleles_coding_sequences(progressive_analysis_instance, alignment_information_instance)
        self.find_matching_combinations_per_position(progressive_analysis_instance, alignment_information_instance)
        self.get_typing_results()

    def get_result_alleles_coding_sequences(self, progressive_analysis_instance, alignment_information_instance):
        """Gets coding sequences of alleles in result_alleles_dictionary attribute from Progressive Analysis instance.
        Sequences are stored in the dictionary, as values for every allele name key.
        Sequences are obtained from the reference file.

        Parameters
        ----------
        progressive_analysis_instance : class ProgressiveAnalysis instance
            Instance inherited from class ProgressiveAnalysis
        alignment_information_instance : class AlignmentInformation instance
            Instance inherited from class AlignmentInformation
        """

        with open(alignment_information_instance.reference.file_name) as alignment_reference:
            for allele in alignment_reference.readlines()[12:-1]:   # First sequence in reference is in line 12
                allele = allele.replace(" ", "")
                kir_gene = allele.split("*")[0]
                if kir_gene == "KIR2DL5A" or kir_gene == "KIR2DL5B":
                    kir_gene = "KIR2DL5"  # KIR genes 2DL5A and 2DL5B are analyzed as a single gene locus

                if kir_gene in progressive_analysis_instance.result_alleles_dictionary.keys():
                    allele_id = allele.split("\t")[0]
                    if allele_id in progressive_analysis_instance.result_alleles_dictionary[kir_gene]:
                        allele_sequence = allele.split("\t")[:-4]  # Last four items in the list are \t characters
                        allele_sequence = allele_sequence[1].replace("|", "")
                        allele_sequence = list(allele_sequence)
                        coding_sequence = []
                        for position in alignment_information_instance.proportion_dictionary.keys():
                            coding_sequence.append(allele_sequence[position])
                        progressive_analysis_instance.result_alleles_dictionary[kir_gene][allele_id] = coding_sequence
        alignment_reference.close()

    def find_matching_combinations_per_position(self, progressive_analysis_instance, alignment_information_instance):
        """Per position in the class AlignmentInfo instance, it finds perfect matching genotype combinations using the
        result alleles combinations coding sequences. Calls methods is_discriminant() and check_genotype_combinations().
        First, it creates all possible genotype combinations of 2 alleles per detected KIR gene.

        Parameters
        ----------
        progressive_analysis_instance : class ProgressiveAnalysis instance
            Instance inherited from class ProgressiveAnalysis
        alignment_information_instance : class AlignmentInformation instance
            Instance inherited from class AlignmentInformation
        """

        result_combinations_dict = {}
        for kir_gene in progressive_analysis_instance.result_alleles_dictionary:
            result_combinations_dict[kir_gene] = (
                list(it.combinations_with_replacement(list(progressive_analysis_instance.result_alleles_dictionary
                                                           [kir_gene].keys()), 2)))
        single_gene_combinations = list(result_combinations_dict.values())
        self.all_genotype_combinations = list(it.product(*single_gene_combinations))

        self.genotype_combinations_per_position = []
        for coding_position, alignment_position in enumerate(alignment_information_instance.proportion_dictionary):
            present_nucleotides = []
            for nucleotide in \
                    alignment_information_instance.proportion_dictionary[alignment_position]["Present Nucleotides"]:
                if nucleotide != ".":
                    present_nucleotides.append(nucleotide)
            if len(present_nucleotides) > 1:
                self.coding_sequences_per_position = {}
                for gene_index, kir_gene in enumerate(progressive_analysis_instance.result_alleles_dictionary):
                    self.coding_sequences_per_position[kir_gene] = {}
                    for allele in progressive_analysis_instance.result_alleles_dictionary[kir_gene]:
                        self.coding_sequences_per_position[kir_gene][allele] = \
                            progressive_analysis_instance.result_alleles_dictionary[kir_gene][allele][coding_position]
                if self.is_discriminant() is True:
                    self.check_genotype_combinations(present_nucleotides)

    def is_discriminant(self):
        """It determines whether a position in the alignment is discriminant. This depends on the number of present
        nucleotides, and the distribution of these nucleotides in the coding sequences of result alleles from
        Progressive Analysis

        """

        for kir_gene in self.coding_sequences_per_position:
            if len(np.unique(list(self.coding_sequences_per_position[kir_gene].values()))) > 1:
                return True

    def check_genotype_combinations(self, present_nucleotides):
        """In a certain discriminant position, it determines which genotype combinations match the alignment information,
        according to the coding_sequences_per_position attribute. It updates the combined_analysis_results with the
        found perfect combinations.

        Parameters
        ----------
        present_nucleotides : list
            List of nucleotides defined as present by Progressive Analysis in a certain position
        """

        combinations_in_position = []
        for combination in self.all_genotype_combinations:
            covered_nucleotides = []
            for single_gene_combination in combination:
                for allele in single_gene_combination:
                    kir_gene = allele.split("*")[0]
                    if kir_gene == "KIR2DL5A" or kir_gene == "KIR2DL5B":
                        kir_gene = "KIR2DL5"  # KIR genes 2DL5A and 2DL5B are analyzed as a single gene locus
                    covered_nucleotides.append(self.coding_sequences_per_position[kir_gene][allele])
            if set(present_nucleotides).issubset(set(covered_nucleotides)) is True:
                combinations_in_position.append(combination)
        self.genotype_combinations_per_position.append(combinations_in_position)

    def get_typing_results(self):
        """Iterates over self.genotype_combinations_per_position to find final typing results;
        genotype combinations that match all discriminant positions in the alignment data.
        Creates a result_alleles_dictionary to compare with allelic results from Progressive analysis

        Raises
        ------
        Combined analysis is not applicable when less than 3 alleles are possible allele results per detected KIR genes

        """

        for final_list_index in range(len(self.genotype_combinations_per_position)):
            if final_list_index == 0:
                self.typing_result = list(set(self.genotype_combinations_per_position[final_list_index]))
            else:
                self.typing_result = set(self.typing_result).intersection(set(self.genotype_combinations_per_position
                                                                              [final_list_index]))
        if self.typing_result is None:
            print("Combined analysis not applicable")

        else:
            remaining_alleles = []
            for combination in self.typing_result:
                for kir_gene in combination:
                    for allele in kir_gene:
                        remaining_alleles.append(allele)
            remaining_alleles = np.unique(remaining_alleles)

            for allele in remaining_alleles:
                kir_gene = allele.split("*")[0]
                if kir_gene == "KIR2DL5A" or kir_gene == "KIR2DL5B":
                    kir_gene = "KIR2DL5"   # KIR genes 2DL5A and 2DL5B are analyzed as a single gene locus
                if kir_gene not in self.result_alleles_dictionary:
                    self.result_alleles_dictionary[kir_gene] = []
                self.result_alleles_dictionary[kir_gene].append(allele)
