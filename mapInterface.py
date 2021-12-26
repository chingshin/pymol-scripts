# -*- coding: utf-8 -*-
"""

This is a python/pymol script to find interface residues between two chains in a complex.

Author: chingshin, 2021/03/04

"""

from pymol import cmd, stored

def mapInterface(sele1, sele2, cutoff=6.0):
    """
    mapInterface -- finds 'interface' residues between two chains in a complex.

    Parameters
    ----------
    sele1 : TYPE: string (pymol selection)
        The first chain in which we search for residues at an interface with 
        the second chain.
    sele2 : TYPE: string (pymol selection)
        The second chain in which we search for residues at an interface with 
        the first chain.
    cutoff : TYPE: float, optional
        The cutoff distance between the first chain and second chain for 
        searching an interface. The default is 6.0.

    Returns
    -------
    None.

    """
    seen = []
    stored.r = []
    interface = cmd.get_unused_name("interface_")
    cmd.select("temp1", f"byresi {sele1} or {sele2} within {cutoff} of {sele1}")
    cmd.select("temp2", f"byresi {sele2} or {sele1} within {cutoff} of {sele2}")
    cmd.select(f"{interface}", "temp1 and temp2")
    cmd.iterate(f"{interface} and name ca", "stored.r.append((model, chain, oneletter, resi))")
    
    for (model, chain, oneletter, resi) in sorted(stored.r, key=lambda tup: int(tup[3])):
        cmd.select("temp3", f"byresi {sele1} or {sele2} within {cutoff} of model {model} and chain {chain} and resi {resi}")
        cmd.select(f"{chain}_{oneletter}{resi}", f"temp3 and (not chain {chain}) or (model {model} and chain {chain} and resi {resi})")
        cmd.group(f"{interface}_{chain}", f"{chain}_{oneletter}{resi}")
    
    cmd.deselect()
    cmd.delete("temp1")
    cmd.delete("temp2")
    cmd.delete("temp3")    
    
cmd.extend("mapInterface", mapInterface)
