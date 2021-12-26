# -*- coding: utf-8 -*-
"""

This is a python/pymol script to assign FR and CDR regions of Fv in an antibody.
This script adapts codes from Xin Yu's annotate_v pymol script (https://pymolwiki.org/index.php/Annotate_v).
The CDR definition (kabat, chothia, imgt, martin and aho) is from https://doi.org/10.3389/fimmu.2018.02278

Author: chingshin, 2021/12/25

"""

import requests
from pymol import cmd

class Fv():
    """
    class `Fv`.

    Initiator `__init__` has 3 parameters:

    :param aaseq: STRING: A SINGLE-LETTER, amino acid sequence corresponding to the complete VH or VL chain.

    :param scheme: STRING: "kabat", "chothia", "martin", "contact", "imgt", or "aho". Must be in LOWERCASE
    
    Class has 2 methods. `retrieve()`: retrieves numbered seqs from Abnum website. `analyze` determine the FR and CDR regions.
    """

    def __init__(self, aaseq, scheme):

        self.aaseq = aaseq
        self.scheme = scheme

    def __repr__(self):
        return "Annotation of VH or VL sequence using Kabat, Chothia, Martin, Contact, IMGT or Aho scheme"

    def analyze(self, abnumPage, scheme):
        """
        Define CDR and FR regions based on the numbered sequence returned from website

        :param chain: STRING, "H" or "L" in uppercase
        :param lst: LIST, a list of residue and their corresponding numbers in kabat or chothia scheme
        :return: DICTIONARY, a dictionary of strings, where the region is the key and the corresponding peptide is the value 
        :raises: `ValueError` if any of the FR or CDR region is missing

        """
        
        lst = abnumPage.text.split()
        chain = lst[0][0]
        
        if chain == "L":
            try:
                if scheme == "kabat":
                    L_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("L24"), 2)])
                    L_CDR1 = "".join([lst[i + 1] for i in range(lst.index("L24"), lst.index("L35"), 2)])
                    L_FR2 = "".join([lst[i + 1] for i in range(lst.index("L35"), lst.index("L50"), 2)])
                    L_CDR2 = "".join([lst[i + 1] for i in range(lst.index("L50"), lst.index("L57"), 2)])
                    L_FR3 = "".join([lst[i + 1] for i in range(lst.index("L57"), lst.index("L89"), 2)])
                    L_CDR3 = "".join([lst[i + 1] for i in range(lst.index("L89"), lst.index("L98"), 2)])
                    L_FR4 = "".join([lst[i + 1] for i in range(lst.index("L98"), len(lst), 2)])

                elif scheme == "chothia":
                    L_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("L26"), 2)])
                    L_CDR1 = "".join([lst[i + 1] for i in range(lst.index("L26"), lst.index("L33"), 2)])
                    L_FR2 = "".join([lst[i + 1] for i in range(lst.index("L33"), lst.index("L50"), 2)])
                    L_CDR2 = "".join([lst[i + 1] for i in range(lst.index("L50"), lst.index("L53"), 2)])
                    L_FR3 = "".join([lst[i + 1] for i in range(lst.index("L53"), lst.index("L91"), 2)])
                    L_CDR3 = "".join([lst[i + 1] for i in range(lst.index("L91"), lst.index("L97"), 2)])
                    L_FR4 = "".join([lst[i + 1] for i in range(lst.index("L97"), len(lst), 2)])
                
                elif scheme == "martin":
                    L_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("L30"), 2)])
                    L_CDR1 = "".join([lst[i + 1] for i in range(lst.index("L30"), lst.index("L36"), 2)])
                    L_FR2 = "".join([lst[i + 1] for i in range(lst.index("L36"), lst.index("L46"), 2)])
                    L_CDR2 = "".join([lst[i + 1] for i in range(lst.index("L46"), lst.index("L56"), 2)])
                    L_FR3 = "".join([lst[i + 1] for i in range(lst.index("L56"), lst.index("L89"), 2)])
                    L_CDR3 = "".join([lst[i + 1] for i in range(lst.index("L89"), lst.index("L97"), 2)])
                    L_FR4 = "".join([lst[i + 1] for i in range(lst.index("L97"), len(lst), 2)])
                                        
                elif scheme == "contact":
                    L_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("L30"), 2)])
                    L_CDR1 = "".join([lst[i + 1] for i in range(lst.index("L30"), lst.index("L37"), 2)])
                    L_FR2 = "".join([lst[i + 1] for i in range(lst.index("L37"), lst.index("L46"), 2)])
                    L_CDR2 = "".join([lst[i + 1] for i in range(lst.index("L46"), lst.index("L56"), 2)])
                    L_FR3 = "".join([lst[i + 1] for i in range(lst.index("L56"), lst.index("L89"), 2)])
                    L_CDR3 = "".join([lst[i + 1] for i in range(lst.index("L89"), lst.index("L97"), 2)])
                    L_FR4 = "".join([lst[i + 1] for i in range(lst.index("L97"), len(lst), 2)])

                elif scheme == "imgt":
                    L_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("L27"), 2)])
                    L_CDR1 = "".join([lst[i + 1] for i in range(lst.index("L27"), lst.index("L39"), 2)])
                    L_FR2 = "".join([lst[i + 1] for i in range(lst.index("L39"), lst.index("L56"), 2)])
                    L_CDR2 = "".join([lst[i + 1] for i in range(lst.index("L56"), lst.index("L66"), 2)])
                    L_FR3 = "".join([lst[i + 1] for i in range(lst.index("L66"), lst.index("L105"), 2)])
                    L_CDR3 = "".join([lst[i + 1] for i in range(lst.index("L105"), lst.index("L118"), 2)])
                    L_FR4 = "".join([lst[i + 1] for i in range(lst.index("L118"), len(lst), 2)])
                
                else:  # aho scheme
                    L_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("L25"), 2)])
                    L_CDR1 = "".join([lst[i + 1] for i in range(lst.index("L25"), lst.index("L41"), 2)])
                    L_FR2 = "".join([lst[i + 1] for i in range(lst.index("L41"), lst.index("L58"), 2)])
                    L_CDR2 = "".join([lst[i + 1] for i in range(lst.index("L58"), lst.index("L78"), 2)])
                    L_FR3 = "".join([lst[i + 1] for i in range(lst.index("L78"), lst.index("L109"), 2)])
                    L_CDR3 = "".join([lst[i + 1] for i in range(lst.index("L109"), lst.index("L138"), 2)])
                    L_FR4 = "".join([lst[i + 1] for i in range(lst.index("L138"), len(lst), 2)])
                
                self.regions = {"L-FR1": L_FR1, "L-CDR1": L_CDR1, "L-FR2": L_FR2, "L-CDR2": L_CDR2, "L-FR3": L_FR3, 
                       "L-CDR3": L_CDR3, "L-FR4":L_FR4}

                return self.regions

            except ValueError:
                print("Unable to define complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured")
        else:
            try:
                if scheme == "kabat":
                    H_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("H31"), 2)])
                    H_CDR1 = "".join([lst[i + 1] for i in range(lst.index("H31"), lst.index("H36"), 2)])
                    H_FR2 = "".join([lst[i + 1] for i in range(lst.index("H36"), lst.index("H50"), 2)])
                    H_CDR2 = "".join([lst[i + 1] for i in range(lst.index("H50"), lst.index("H66"), 2)])
                    H_FR3 = "".join([lst[i + 1] for i in range(lst.index("H66"), lst.index("H95"), 2)])
                    H_CDR3 = "".join([lst[i + 1] for i in range(lst.index("H95"), lst.index("H103"), 2)])
                    H_FR4 = "".join([lst[i + 1] for i in range(lst.index("H103"), len(lst), 2)])

                elif scheme == "chothia":
                    H_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("H26"), 2)])
                    H_CDR1 = "".join([lst[i + 1] for i in range(lst.index("H26"), lst.index("H33"), 2)])
                    H_FR2 = "".join([lst[i + 1] for i in range(lst.index("H33"), lst.index("H52"), 2)])
                    H_CDR2 = "".join([lst[i + 1] for i in range(lst.index("H52"), lst.index("H57"), 2)])
                    H_FR3 = "".join([lst[i + 1] for i in range(lst.index("H57"), lst.index("H95"), 2)])
                    H_CDR3 = "".join([lst[i + 1] for i in range(lst.index("H95"), lst.index("H103"), 2)])
                    H_FR4 = "".join([lst[i + 1] for i in range(lst.index("H103"), len(lst), 2)])

                elif scheme == "martin":
                    H_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("H30"), 2)])
                    H_CDR1 = "".join([lst[i + 1] for i in range(lst.index("H30"), lst.index("H36"), 2)])
                    H_FR2 = "".join([lst[i + 1] for i in range(lst.index("H36"), lst.index("H47"), 2)])
                    H_CDR2 = "".join([lst[i + 1] for i in range(lst.index("H47"), lst.index("H59"), 2)])
                    H_FR3 = "".join([lst[i + 1] for i in range(lst.index("H59"), lst.index("H95"), 2)])
                    H_CDR3 = "".join([lst[i + 1] for i in range(lst.index("H95"), lst.index("H102"), 2)])
                    H_FR4 = "".join([lst[i + 1] for i in range(lst.index("H102"), len(lst), 2)])
                
                elif scheme == "contact":
                    H_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("H30"), 2)])
                    H_CDR1 = "".join([lst[i + 1] for i in range(lst.index("H30"), lst.index("H36"), 2)])
                    H_FR2 = "".join([lst[i + 1] for i in range(lst.index("H36"), lst.index("H47"), 2)])
                    H_CDR2 = "".join([lst[i + 1] for i in range(lst.index("H47"), lst.index("H59"), 2)])
                    H_FR3 = "".join([lst[i + 1] for i in range(lst.index("H59"), lst.index("H93"), 2)])
                    H_CDR3 = "".join([lst[i + 1] for i in range(lst.index("H93"), lst.index("H102"), 2)])
                    H_FR4 = "".join([lst[i + 1] for i in range(lst.index("H102"), len(lst), 2)])

                elif scheme == "imgt":
                    H_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("H27"), 2)])
                    H_CDR1 = "".join([lst[i + 1] for i in range(lst.index("H27"), lst.index("H39"), 2)])
                    H_FR2 = "".join([lst[i + 1] for i in range(lst.index("H39"), lst.index("H56"), 2)])
                    H_CDR2 = "".join([lst[i + 1] for i in range(lst.index("H56"), lst.index("H66"), 2)])
                    H_FR3 = "".join([lst[i + 1] for i in range(lst.index("H66"), lst.index("H105"), 2)])
                    H_CDR3 = "".join([lst[i + 1] for i in range(lst.index("H105"), lst.index("H118"), 2)])
                    H_FR4 = "".join([lst[i + 1] for i in range(lst.index("H118"), len(lst), 2)])
                
                else:  # aho scheme
                    H_FR1 = "".join([lst[i + 1] for i in range(0, lst.index("H25"), 2)])
                    H_CDR1 = "".join([lst[i + 1] for i in range(lst.index("H25"), lst.index("H41"), 2)])
                    H_FR2 = "".join([lst[i + 1] for i in range(lst.index("H41"), lst.index("H58"), 2)])
                    H_CDR2 = "".join([lst[i + 1] for i in range(lst.index("H58"), lst.index("H78"), 2)])
                    H_FR3 = "".join([lst[i + 1] for i in range(lst.index("H78"), lst.index("H109"), 2)])
                    H_CDR3 = "".join([lst[i + 1] for i in range(lst.index("H109"), lst.index("H138"), 2)])
                    H_FR4 = "".join([lst[i + 1] for i in range(lst.index("H138"), len(lst), 2)])

                self.regions = {"H-FR1": H_FR1, "H-CDR1": H_CDR1, "H-FR2": H_FR2, "H-CDR2": H_CDR2, "H-FR3": H_FR3, 
                       "H-CDR3": H_CDR3, "H-FR4":H_FR4}
                return self.regions

            except ValueError:
                print("Unable to define complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occured in the `analyze()` method")

    def retrieve(self):
        """
        Retrieve numbered residues from Abnum website, and assign the webpage (text) to self.abnum

        :raises: `ValueError` if input scheme is not among "kabat", "chothia", "contact", and "imgt"

        """

        abnum_url = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi"
        scheme_dict = {"chothia":"-c", "kabat":"-k", "imgt":"-imgt", "contact":"-c", "martin":"-m", "aho":"-aho"}
        
        if self.scheme not in ["kabat", "chothia", "contact", "imgt", "martin", "aho"]:
            raise Exception("Incorrect scheme mode. Must be one of the following (lowercase): kabat, chothia, contact, imgt, martin, aho")
             
        query = {"plain": 1, "scheme": scheme_dict.get(self.scheme), "aaseq": self.aaseq}
        myPage = requests.get(abnum_url, params=query)
            
        if len(myPage.text.split()) < 1:
            raise Exception("An error occured in the `retrieve()` method. Check input sequence")
        else:
            self.abnum = myPage
            
def annotateFv(selection, scheme):
    aaseq="".join(cmd.get_fastastr(selection).split("\n")[1:])
    obj=Fv(aaseq, scheme)
    obj.retrieve()
    obj.analyze(obj.abnum, obj.scheme)
    for i in obj.regions.keys():
        cmd.select(i, f"pepseq {obj.regions[i]}")
    cmd.deselect()
    print(f"Annotation in {obj.scheme} scheme:")
    for region, aa in obj.regions.items():
        print(f"{region}: {aa}")

cmd.extend("annotateFv", annotateFv)
