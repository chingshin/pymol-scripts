"""

This is a python/pymol script to analyze a interaction interface between two 
binding partners in great detail.

Given two selections defining an interface, this script finds salt bridges, 
hydrogen bonds, water molecules mediating the interaction, and calculates the 
BSA of each residue in the interface due to each residue on the other side.
For spead, BSA calculation only considers residues within 5 A between two 
binding partners.

Outputs are 1) a csv file (default name: interface_analysis.csv) listing a 2D 
interaction map with BSA annotations, and all electrostatic interactions. 2)
pymol objects that renders BSA calculation (default: percentage of buried) on 
surface and electrostatic interaction network.

Author: chingshin, 2023/02/24
"""
from pymol import cmd
from pymol import util
from pymol import stored

class ResidueInfo:
    def __init__(self, obj, chain, resi):
        self.obj = obj
        self.chain = chain
        self.resi = resi
        stored.resn = []
        cmd.iterate(f"{self.obj} and chain {self.chain} and resi {self.resi} and name CA", "stored.resn.append([resn, oneletter])")
        self.resn = stored.resn[0][0]
        self.oneletter = stored.resn[0][1]
        del stored.resn      
        self.saltBridgePartners=[]
        self.hbondPartners=[]
        self.sharedWaters=[]  #each entry will have ch,res of partner aa sharing water followed by ch,res of water
        self.bsa=0
        self.pctBuried=0
        self.bsaPartners={}

def findSaltBridges(ligand, receptor):
    print(f"find salt bridges between {ligand} and {receptor}")
    cmd.select("bkbn", f"({ligand} or {receptor}) and name CA+C+O+N")
    cmd.select("acid1",f"{ligand} and resn GLU+ASP and elem O and not bkbn")
    cmd.select("acid2",f"{receptor} and resn GLU+ASP and elem O and not bkbn")
    cmd.select("base1",f"{ligand} and resn ARG+LYS+HIS and elem N and not bkbn")
    cmd.select("base2",f"{receptor} and resn ARG+LYS+HIS and elem N and not bkbn")
    cmd.select("salt1","(byres acid1 within 4 of base2) and name CA")
    cmd.select("salt2","(byres base1 within 4 of acid2) and name CA")
    stored.saltAtoms=[]
    cmd.iterate("salt1","stored.saltAtoms.append([chain,resi])")
    cmd.iterate("salt2","stored.saltAtoms.append([chain,resi])")
    saltBridges = []
    for chain, resi in stored.saltAtoms:
        print(chain, resi, end=' ')
        stored.partner=[]
        for p1, p2 in [('base2', 'acid1'),('acid2', 'base1')]:
        #get partner
        #FIX THIS Need to sep Acid and Base!!
            cmd.select("temp", f"(byres {p1} within 4 of (chain {chain} and resi {resi} and {p2})) and name CA")
            cmd.iterate("temp", "stored.partner.append([chain,resi])")
        for pchn,pres in stored.partner:
            if [chain,resi,pchn,pres] in saltBridges:
                continue
            print(pchn,pres)
            #store in residueInfo
            saltBridges.append([chain,resi,pchn,pres])
    
    cmd.delete("bkbn acid1 acid2 base1 base2 salt1 salt2 temp")
    #del stored.saltAtoms
    #del stored.partner
    
    return saltBridges
    
def findHydrogenBonds(ligand, receptor):
    
    print(f"find hydrogen bonds between {ligand} and {receptor}")
    cmd.h_add(f"{ligand} or {receptor}")
    cmd.select("don", f"({ligand} or {receptor}) and (elem n+o and (neighbor hydro))")
    cmd.select("acc", f"({ligand} or {receptor}) and (elem o or (elem n and not (neighbor hydro)))")
    cmd.select("hbond1", f"({ligand} and don) within 3.5 of ({receptor} and acc)")
    cmd.select("hbond2", f"({ligand} and acc) within 3.5 of ({receptor} and don)")
    stored.hbondAtoms = []
    cmd.iterate("hbond1 or hbond2", "stored.hbondAtoms.append([model,chain,resi,name])")
    hbonds = []
    for model, chain, resi, name in stored.hbondAtoms:
        print(chain, resi, name)
        cmd.select("temp", f"({receptor} and (don or acc)) within 3.5 of ({model} and chain {chain} and resi {resi} and name {name})")
        stored.partner=[]
        cmd.iterate("temp", "stored.partner.append([chain, resi, name])")
        for pchn, pres, pname in stored.partner:
            print("  ", pchn, pres, pname)
            if [chain,resi,pchn,pres] in hbonds:
                continue  #multiple interactions between a single residue pair
            hbonds.append([chain, resi, pchn, pres])
            
    
    cmd.delete("don acc hbond1 hbond2 temp")
    cmd.remove("hydro")
    #del stored.hbondAtoms
    #del stored.partner
    
    return hbonds
    
def findInterfaceWaters(ligand, receptor, waterChain="all"):
    print(f"find water molecules mediating interaction between {ligand} and {receptor}")
    
    cmd.h_add(f"{ligand} or {receptor}")
    cmd.select("don", f"({ligand} or {receptor}) and (elem n+o and (neighbor hydro))")
    cmd.select("acc", f"({ligand} or {receptor}) and (elem o or (elem n and not (neighbor hydro)))")
    cmd.select("sharedWat", f"(solvent within 3.5 of ({ligand} and (don or acc))) and (solvent within 3.5 of ({ligand} and (don or acc))) and {waterChain}")
    stored.sharedWaters = []
    cmd.iterate("sharedWat", "stored.sharedWaters.append([model, chain, resi])")
    sharedWaters=[]
    for wm, wc, wr in stored.sharedWaters:
        cmd.select("temp1", f"(byres ({ligand} and (don or acc)) within 3.5 of {wm} and chain {wc} and resi {wr}) and name CA")
        stored.partner1=[]
        cmd.iterate("temp1", "stored.partner1.append([chain,resi])")
        cmd.select("temp2", f"(byres ({receptor} and (don or acc)) within 3.5 of {wm} and chain {wc} and resi {wr}) and name CA")
        stored.partner2=[]
        cmd.iterate("temp2", "stored.partner2.append([chain,resi])")
        for chain1, resi1 in stored.partner1:
            for chain2, resi2 in stored.partner2:
                if [chain1, resi1, chain2, resi2, wc, wr] in sharedWaters:
                    continue #avoid duplicates
                sharedWaters.append([chain1, resi1, chain2, resi2, wc, wr])
    
    cmd.delete("don acc temp1 temp2 sharedWat")
    cmd.remove("hydro")
    #del stored.sharedWaters
    #del stored.partner1
    #del stored.partner2
    print(sharedWaters)
    
    return sharedWaters

def calculateBSA(wholeComplex, ligand, receptor, cutoff=5.0, includeSharedWaters=True, output="BSA_report.txt"):
    temp = cmd.get_unused_name("temp_") #use this to get rid of pymol memory problem
    #set dot_solvent to calculate SASA and use higer sampling density
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", 3)
    
    #select residues to be included in calculation
    stored.interDict ={}
    cmd.iterate(f"(byres {wholeComplex} and {ligand} within {cutoff} of {wholeComplex} and {receptor}) and name CA", f"stored.interDict[(chain, resi)] = []")
    for ligSel in stored.interDict:
        cmd.iterate(f"(byres {wholeComplex} and {receptor} within {cutoff} of {wholeComplex} and {ligand} and chain {ligSel[0]} and resi {ligSel[1]}) and name CA", f"stored.interDict[{ligSel}].append([chain, resi])")
    
    print(f"calculate BSA for {[i[0]+i[1] for i in stored.interDict.keys()]}")
    
    #make temporary models depending on include shared water or not
    if includeSharedWaters:
        cmd.h_add(f"{wholeComplex} and ({ligand} or {receptor})")
        cmd.select("don", f"({wholeComplex} and ({ligand} or {receptor})) and (elem n+o and (neighbor hydro))")
        cmd.select("acc", f"({wholeComplex} and ({ligand} or {receptor})) and (elem o or (elem n and not (neighbor hydro)))")
        cmd.select("sharedW", f"(({wholeComplex} and solvent) within 3.5 of ({ligand} and (don or acc))) and (({wholeComplex} and solvent) within 3.5 of ({receptor} and (don or acc)))")
        cmd.remove("hydro")
        
        print("include shared water:", end = " ")
        cmd.iterate("sharedW", "print(chain+resi, end = ' ')")
        print("")
        cmd.select("tempLig", f"{wholeComplex} and ({ligand} or sharedW)")
        cmd.select("tempCplx", f"{wholeComplex} and ({receptor} or {ligand} or sharedW)")
        cmd.flag('ignore','all','clear') #make sure waters aren't ignored
    
    else:
        cmd.select("tempLig", f"{wholeComplex} and {ligand}")
        cmd.select("tempCplx", f"{wholeComplex} and ({receptor} or {ligand})")
        
    BSAlist = []
    #BSAreport = []
    out = open(output, "w")
    print("chn resi uASA bASA BSA %b BSAPrtnr scale Partners", file = out)    
    for ligSel in stored.interDict:
        cmd.flag("ignore", "not tempLig")
        unboundASA = cmd.get_area(f"tempLig and chain {ligSel[0]} and resi {ligSel[1]}")
        cmd.flag("ignore", "not tempCplx")
        boundASA = cmd.get_area(f"tempCplx and chain {ligSel[0]} and resi {ligSel[1]}")
        BSA = unboundASA - boundASA
        pctBuried = 100*BSA/unboundASA
        totalRecBSA = 0
        recBsaList = []
        for recRes in stored.interDict[ligSel]:
            cmd.flag("ignore", "tempCplx", "clear")
            cmd.select(f"{temp}", f"({wholeComplex} and {receptor} and chain {recRes[0]} and resi {recRes[1]}) or tempLig")
            cmd.flag("ignore", f"not {temp}")
            partialBoundASA = cmd.get_area(f"{temp} and chain {ligSel[0]} and resi {ligSel[1]}")
            partialBSA = unboundASA - partialBoundASA
            recBsaList.append([recRes[0], recRes[1], partialBSA])
            totalRecBSA += partialBSA
        if totalRecBSA > 0.0:
            scaleFactor = BSA/totalRecBSA
            for i in recBsaList:
                i[2] = f"{i[2]*scaleFactor:.2f}"
        else:
            scaleFactor = 1
        
        ligBsaList = [ligSel[0], ligSel[1], f"{unboundASA:.2f}", f"{boundASA:.2f}", f"{BSA:.2f}", f"{pctBuried:.2f}", f"{totalRecBSA:.2f}", f"{scaleFactor:.2f}"]
        BSAlist.append(ligBsaList + recBsaList)
        #flatrecBsaList = [item for sublist in recBsaList for item in sublist] #flaten (remove sublists in) recBsaList
        #BSAreport.append(ligBsaList + flatrecBsaList)
    
        #write an BSA report
        print(*ligBsaList, end = " ", file = out)
        for i in recBsaList:
            print(*i, end = " ", file = out)
        print("", file = out)
    out.close() 
    
    print(f"Please see results in {output}")
    
    #clean up
    cmd.flag("ignore", "all", "clear")
    cmd.delete("don acc sharedW")
    #cmd.delete("temp*")
    #del stored.interDict
    
    return BSAlist

def findCloseContacts(ligand, receptor, contactDistance=3.8):
    cmd.select("tempClose1", f"byresi {ligand} within {contactDistance} of {receptor}")
    stored.contactDict = {}
    cmd.iterate("tempClose1 and name ca", "stored.contactDict.update({(chain, resi):[]})")
    for k in stored.contactDict:
        cmd.select("tempClose2", f"byresi {receptor} within {contactDistance} of {ligand} and chain {k[0]} and resi {k[1]}")
        cmd.iterate("tempClose2 and name ca", f"stored.contactDict[{k}].append([chain, resi])")
    
    contactList = [[*k, *stored.contactDict[k][v]] for k in stored.contactDict for v in range(len(stored.contactDict[k]))]
    
    cmd.delete("tempClose*")
    #del stored.contactDict
    
    print(f"find residue pairs between {ligand} and {receptor} within {contactDistance} Angstrom of contact distance:")
    print(*contactList, sep = "\n")
    
    return contactList
        
def analyzeInterface(wholeComplexObj, ligand, receptor, cutoff=5.0, contactDistance=3.8, includeSharedWaters=True, output="interface_analysis.csv", usePctBuried=True):
    
    '''
    Description
    ---------------------------------------------------------------------------
    "analyzeInterface" finds all interactions between two specified selections 
    in a pymol object.  This includes salt bridges, hydrogen bonds, shared 
    waters, and the pairwise buried surface for each pair of residues that 
    interact at the interface. The pairwise buried surface is calculated by 
    determining the solvent accessible surface area that is buried in one 
    binding partner by each residue individually on the other. These areas
    are then normalized so that the sum of the areas buried on each residue by 
    individual interacting amino acids is the same as the area buried on that
    residue by the formation of the entire complex.
    
    Parameters
    ---------------------------------------------------------------------------
    wholeComplexObj: The original structure object to be analyzed
    sel1: One binding partner in pymol selection string
    sel2: The other binding partner in pymol selection string
    cutoff: Residues in one partner within cutoff from the other partner will 
            be considered in each pairwise calculation
    contactDistance: Distance between atoms to be considered as being in 
                     contact
    includeSharedWaters: If True, the area buried by a water coordinated by 
                         both binding partners will be included in 
                         calculation of the pairwise BSA of each coordinating 
                         residue
    output: the name of the csv file to be generated after analysis
    usePctBuried: if True, surface is colored based on percentage of buried 
                  in the complex. otherwise, surface is colored based on BSA
                  
    Outputs
    ---------------------------------------------------------------------------
    1) a csv file (default name: interface_analysis.csv) in the current 
    directory listing a 2D interaction map with BSA annotations, and all 
    electrostatic interactions
    
    2) Pymol objects are created including:
        Each binding partner as a separate object. These are displayed as 
        surfaces colored by per-residue area buried by complex formation. White 
        residues were not considered, green have least BSA and red have most 
        BSA.
        
        interface water molecules mediating the interaction are shown as 
        licorice.
        
        Hydrogen bonds and salt bridges are shown as distance objects.

    '''
    
    print(wholeComplexObj, ligand, receptor, cutoff, includeSharedWaters)
    
    #make selections for the analysis
    sel1 = f"(not hetatm and not solvent and {ligand} and {wholeComplexObj})"
    sel2 = f"(not hetatm and not solvent and {receptor} and {wholeComplexObj})"
    quickWat = f"((solvent within 3.5 of {sel1}) and (solvent within 3.5 of {sel2}) and {wholeComplexObj})"
    
    #select the residues that will be included in calculation
    cmd.select("nearRes1", f" byres {sel1} within {cutoff} of {sel2}")
    cmd.select("nearRes2", f"byres {sel2} within {cutoff} of {sel1}")
    cmd.select("nearResCplx", f"byres nearRes1 or nearRes2")
    
    #get the selected residues into python
    stored.residues1 = []
    stored.residues2 = []
    cmd.iterate("nearRes1 and name ca", "stored.residues1.append([model, chain, resi])")
    cmd.iterate("nearRes2 and name ca", "stored.residues2.append([model, chain, resi])")
 
    #create objects to store residue properties
    residueInfo1 = {}
    for (model, chain, resi) in stored.residues1:
        residueInfo1[chain+resi] = ResidueInfo(model, chain, resi)

    residueInfo2 = {}
    for (model, chain, resi) in stored.residues2:
        residueInfo2[chain+resi] = ResidueInfo(model, chain ,resi)
    
    #find salt bridges in the interaction interface and load salt bridge pairs to residueInfo
    sbPairs = findSaltBridges(sel1, sel2)
    for i in sbPairs:
        residueInfo1[i[0]+i[1]].saltBridgePartners.append([i[2], i[3]])
        residueInfo2[i[2]+i[3]].saltBridgePartners.append([i[0], i[1]])
        
    #find hydrogen bonds in the interaction interface and load hydrogen-bond pairs to residueInfo
    hbPairs = findHydrogenBonds(sel1, sel2)
    for i in hbPairs:
        residueInfo1[i[0]+i[1]].hbondPartners.append([i[2], i[3]])
        residueInfo2[i[2]+i[3]].hbondPartners.append([i[0], i[1]])
        
    #find shared water moleucules in the interaction interface and load them to residueInfo
    intfWaters = findInterfaceWaters(sel1, sel2, quickWat)
    for i in intfWaters:
        #pairs sharing water might be outside distance cutoff so add them
        if i[0]+i[1] not in residueInfo1:
            residueInfo1[i[0]+i[1]] = ResidueInfo(wholeComplexObj, i[0], i[1])
            stored.residues1.append([wholeComplexObj, i[0], i[1]])
        if i[2]+i[3] not in residueInfo2:
            residueInfo2[i[2]+i[3]] = ResidueInfo(wholeComplexObj, i[2], i[3])
            stored.residues2.append([wholeComplexObj, i[2], i[3]])
        residueInfo1[i[0]+i[1]].sharedWaters.append([i[2], i[3], i[4], i[5]])
        residueInfo2[i[2]+i[3]].sharedWaters.append([i[0], i[1], i[4], i[5]])   

    #calculate BSA of each binding partner in the interface
    sel1BsaList = calculateBSA(wholeComplexObj, sel1, sel2, cutoff, includeSharedWaters, "sel1_BSA_report.txt")
    sel2BsaList = calculateBSA(wholeComplexObj, sel2, sel1, cutoff, includeSharedWaters, "sel2_BSA_report.txt")
    for i in sel1BsaList:
        residueInfo1[i[0]+i[1]].bsa = i[4]
        residueInfo1[i[0]+i[1]].pctBuried = i[5]
        residueInfo1[i[0]+i[1]].bsaPartners = {a[0]:a[1] for a in [(b[0]+b[1], b[2]) for b in i[8:]]}
    for i in sel2BsaList:
        residueInfo2[i[0]+i[1]].bsa = i[4]
        residueInfo2[i[0]+i[1]].pctBuried = i[5]
        residueInfo2[i[0]+i[1]].bsaPartners = {a[0]:a[1] for a in [(b[0]+b[1], b[2]) for b in i[8:]]}
    
    #generate close residue pairs using contactDistance defined by input (default: 3.8 Angstrom)
    ccList = findCloseContacts(sel1, sel2, contactDistance)
    
    #generate and output results to a csv file
    csv = open(output, "w")
    nameLine = " , ,"
    totalBSALine = " ,total,"
    for k2 in residueInfo2:
        nameLine += f"{residueInfo2[k2].chain} {residueInfo2[k2].oneletter}{residueInfo2[k2].resi},"
        totalBSALine += f"{residueInfo2[k2].bsa},"
    print(nameLine, "\n", totalBSALine, file=csv)
    
    sumPairwiseBSA = 0.0
    for k1 in residueInfo1:
        print(f"{residueInfo1[k1].chain} {residueInfo1[k1].oneletter}{residueInfo1[k1].resi}, {residueInfo1[k1].bsa},", end = " ", file = csv)
        for k2 in residueInfo2:
            pairwiseBSA = 0.0
            electroStatic = ""
            #calculate pairwiseBSA
            if k2 in residueInfo1[k1].bsaPartners:
                pairwiseBSA += float(residueInfo1[k1].bsaPartners[k2])
            if k1 in residueInfo2[k2].bsaPartners:
                pairwiseBSA += float(residueInfo2[k2].bsaPartners[k1])
            #annotate electroStatic interaction
            if k2 in [i[0]+i[1] for i in residueInfo1[k1].saltBridgePartners]:
                electroStatic += "S"
            if k2 in [i[0]+i[1] for i in residueInfo1[k1].hbondPartners]:
                electroStatic += "H"
            if k2 in [i[0]+i[1] for i in residueInfo1[k1].sharedWaters]:
                electroStatic += "W"
            
            #print pairwiseBSA and electroStatic info to the csv file
            if pairwiseBSA > 0.0:
                print(f"{pairwiseBSA:.2f} {electroStatic},", end = " ", file = csv)
                sumPairwiseBSA += pairwiseBSA
            else:
                print(f" {electroStatic},", end = " ", file = csv)
        print("", file = csv)
    print("", file = csv)
    print(f"Sum of pairwise BSA,{sumPairwiseBSA:.2f}", file = csv)
    cmd.flag("ignore", f"not {sel1}")
    totalAreaFree1 = cmd.get_area(f"{sel1}")
    cmd.flag("ignore", f"not {sel2}")
    totalAreaFree2 = cmd.get_area(f"{sel2}")
    cmd.flag("ignore", f"not ({sel1} or {sel2} or {quickWat})")
    totalAreaCplx = cmd.get_area(f"{sel1} or {sel2} or {quickWat}")
    cmd.flag("ignore", "all", "clear")
    interfaceBSA = totalAreaFree1 + totalAreaFree2 - totalAreaCplx
    print(f"Full interface BSA,{interfaceBSA}", file = csv)
    ##summary of salt bridges
    print("Salt Bridged Residue Pairs", file = csv)
    print("Chain,Resi,Chain,Resi", file = csv)
    for s in sbPairs:
        print(*s, sep = ",", file = csv)
    ##summary of hydrogen bonds
    print("Hydrogen Bonded Residue Pairs", file = csv)
    print("Chain,Resi,Chain,Resi", file = csv)
    for h in hbPairs:
        print(*h, sep = ",", file = csv)
    ##summary of shared water molecules
    print("Shared Water Residue Pairs", file = csv)
    print("Chain,Resi,Chain,Resi,solChn,solResi", file = csv)
    for w in intfWaters:
        print(*w, sep = ",", file = csv)
    ##filter "core residues" using BSA > 20 or electrostatic interaction
    print("Residues with BSA > 20 or electrostatics", file = csv)
    print("Chain,Resi,Type,BSA,%Buried,Electrostatic Interactions", file = csv)
    allResiInfo = {**residueInfo1, **residueInfo2}
    for k in allResiInfo:
        if float(allResiInfo[k].bsa) > 20.0 or allResiInfo[k].hbondPartners or allResiInfo[k].saltBridgePartners or allResiInfo[k].sharedWaters:
            print(f"{allResiInfo[k].chain},{allResiInfo[k].resi},{allResiInfo[k].resn},{allResiInfo[k].bsa},{allResiInfo[k].pctBuried}", end = " ", file = csv)
        else:
            continue
            
        if allResiInfo[k].hbondPartners:
                print(",H", end=' ', file=csv)
        if allResiInfo[k].saltBridgePartners:
                print(",S", end=' ', file=csv)
        if allResiInfo[k].sharedWaters:
                print(",W", end=' ', file=csv)
        print(file=csv)
    ##summary of residue pairs within defined contact distance   
    print(f"Residues having atoms within {contactDistance} Angstrom", file = csv)
    print("chain,resi,chain,resi", file = csv)
    for i in ccList:
        print(*i, sep = ",", file = csv)
    csv.close()
    
    #display properties in pymol
    ##create new objects for visuallization
    intObj1 = cmd.get_unused_name("intLig")
    intObj2 = cmd.get_unused_name("intRec")
    intCplx = cmd.get_unused_name("intCplx")
    cmd.create(f"{intCplx}", f"{sel1} or {sel2} or {quickWat}")
    cmd.create(f"{intObj1}", f"{sel1}")
    cmd.create(f"{intObj2}", f"{sel2}")
    cmd.alter(f"{intObj1} or {intObj2}", "q=0.0")  #set occs to 0, will store BSA
    cmd.flag("ignore", "all", "clear") #make sure waters aren't ignored
    updatedSel1 = []
    updatedSel2 = []
    for i in residueInfo1:
       updatedSel1.append(f"chain {residueInfo1[i].chain} and resi {residueInfo1[i].resi}")
    for i in residueInfo2:
       updatedSel2.append(f"chain {residueInfo2[i].chain} and resi {residueInfo2[i].resi}")
    cmd.select(f"{intObj1}Face", f"{intObj1} and ({' or '.join(updatedSel1)})")
    cmd.select(f"{intObj2}Face", f"{intObj2} and ({' or '.join(updatedSel2)})")
    cmd.select(f"{intCplx}Face", f"{intCplx} and ({' or '.join(updatedSel1)} or {' or '.join(updatedSel2)})")
    
    bondList = []
    ##make sb selections on intCplx
    cmd.set("h_bond_cutoff_center", "4.0") #increase cutoff for making sb and hb object in pymol
    cmd.set("h_bond_cutoff_edge", "4.0") #increase cutoff for making sb and hb object in pymol
    cmd.set("h_bond_max_angle", "90.0") #increase cutoff for making sb and hb object in pymol
    for i in sbPairs:
        bond = "".join(i)
        cmd.distance(f"sb{bond}", f"{intCplx} and chain {i[0]} and resi {i[1]}", f"{intCplx} and chain {i[2]} and resi {i[3]}", cutoff="4.0", mode="2")
        bondList.append(f"sb{bond}")    
    ##make hb selections on intCplx
    for i in hbPairs:
        bond = "".join(i)
        cmd.distance(f"hb{bond}", f"{intCplx} and chain {i[0]} and resi {i[1]}", f"{intCplx} and chain {i[2]} and resi {i[3]}", cutoff="3.5", mode="2")  
        bondList.append(f"hb{bond}")
    ##make interface water selections on intCplx
    uniqueWater = []
    wNo = 0
    for i in intfWaters:
        if i[4]+i[5] not in uniqueWater:
            uniqueWater.append(i[4]+i[5])
            wNo +=1
            cmd.select(f"w{wNo}{i[4]}{i[5]}", f"{intCplx} and chain {i[4]} and resi {i[5]}")
        cmd.distance(f"w{wNo}{i[0]}{i[1]}{i[2]}{i[3]}", f"{intCplx} and chain {i[4]} and resi {i[5]}", f"{intCplx} and (chain {i[0]} and resi {i[1]} or chain {i[2]} and resi {i[3]})", cutoff="3.5", mode="2")
        bondList.append(f"w{wNo}{i[0]}{i[1]}{i[2]}{i[3]}")
        
    cmd.set("h_bond_cutoff_center", "3.6") #set it back to default
    cmd.set("h_bond_cutoff_edge", "3.2") #set it back to default
    cmd.set("h_bond_max_angle", "63.0") #set it back to default
    
    ##load BSA or percentage of buried ASA to intObj1 or intObj2 for surface visuallization
    if usePctBuried:
        for k1 in residueInfo1:
            cmd.alter(f"{intObj1} and chain {residueInfo1[k1].chain} and resi {residueInfo1[k1].resi}", f"q = {residueInfo1[k1].pctBuried}")
        for k2 in residueInfo2:
            cmd.alter(f"{intObj2} and chain {residueInfo2[k2].chain} and resi {residueInfo2[k2].resi}", f"q = {residueInfo2[k2].pctBuried}")
    else:
        for k1 in residueInfo1:
            cmd.alter(f"{intObj1} and chain {residueInfo1[k1].chain} and resi {residueInfo1[k1].resi}", f"q = {residueInfo1[k1].bsa}")
        for k2 in residueInfo2:
            cmd.alter(f"{intObj2} and chain {residueInfo2[k2].chain} and resi {residueInfo2[k2].resi}", f"q = {residueInfo2[k2].bsa}")
    cmd.hide("everything")
    cmd.util.cbc(f"{intCplx}")
    cmd.util.cnc(f"{intCplx}")
    cmd.show("cartoon", f"{intCplx}")
    cmd.show("sticks", f"{intCplx}Face")
    cmd.show("licorice", "w*")
    for i in bondList:
        cmd.show("dash", i)
        if i[0] == "s": 
            cmd.color("blue", i)
        else: 
            cmd.color("yellow", i)
    cmd.show("surface", f"{intObj1} or {intObj2}")
    cmd.color("white", f"{intObj1} or {intObj2}")
    cmd.spectrum("q", "green_yellow_red", f"byres {intObj1}Face")
    cmd.spectrum("q", "green_yellow_red", f"byres {intObj2}Face")
    
    cmd.delete("temp*")
    
#load analyzeInterface script to pymol env
cmd.extend("analyzeInterface", analyzeInterface)