import sys


def parse_list(text, ch): 
    return map(str.strip, text.split(ch))

def clean(text):
    return text.replace('-', '')


def parser_kegg(file_path, annotations):
    print('Loading kegg annotation...', file=sys.stderr)

    with open(file_path, "r") as file:
        for line in file:
        
            #(seqid, knumber, gsymbol, gname, brite, cazy, cog, disease, drug, enzyme, go, 
            #    module, network, pathway, pubmed, rclass, reaction, tc, desc) = map(str.strip, line.rstrip('\n').split('\t'))


    
            
            (seqid, kegg_seq, knumber, pathway, lgsymbol, lgname, lbrite, lenzyme, lmodule) = map(str.strip, line.rstrip('\n').split('\t'))
    
            if clean(knumber):
                annotations['knumber'][seqid] = list(parse_list(knumber, ','))   
                #seq2knumber[seqid] = frozenset(parse_list(knumber, ','))    
    
            # if clean(reaction):
                # annotations['reaction'][seqid] =list(parse_list(reaction, ','))    
            
            if clean(pathway):
                annotations['pathway'][seqid] = list(parse_list(pathway, ','))
    
            if clean(lmodule): 
                annotations['module'][seqid] = list(parse_list(lmodule, ','))
    
            if clean(lenzyme):
                annotations['enzyme'][seqid] =  list(parse_list(lenzyme, ','))

            if clean(lbrite): 
                annotations['brite'][seqid] = list(parse_list(lbrite, ','))
    
            if clean(lgsymbol):        
                gsymbol_set = set()
                for element in parse_list(lgsymbol, ','):
                    gsymbol_set.update(parse_list(element, ','))    
                annotations['gsymbol'][seqid] = list(gsymbol_set)
    
            if clean(lgname): 
                annotations['gname'][seqid] = list(parse_list(lgname, ','))
    
            # if clean(desc): 
                # annotations['desc'][seqid] = list(parse_list(desc, '|'))

    
    return annotations



def parser_best_desc(file_path, annotations):
    with open(file_path, "r") as file:
        for line in file:
        
            (seqid, descrip, source) = map(str.strip, line.rstrip('\n').split('\t'))
            annotations['best_description'][seqid] = [descrip]
            
    return annotations


def parser_best_name(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:
        
            (seqid, bname, source) = map(str.strip, line.rstrip('\n').split('\t'))
            annotations['bestname'][seqid] = [bname.upper()]
    
    return annotations


def parser_bigg(file_path, annotations):
    with open(file_path, "r") as file:
        for line in file:
        
            seqid, bigg_geneid, bigg_genename, bigg_org, bigg_reacid, bigg_reacname  = map(str.strip, line.rstrip('\n').split('\t'))
            if clean(bigg_reacname): 
                annotations['bigg'][seqid] = bigg_reacname.split(',')
    
    return annotations

    

def parser_card(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:
    
            #seqid, cardname, card_prot, card_ar, card_species = map(str.strip, line.split('\t'))
            seqid, cardname = map(str.strip, line.split('\t'))
            if clean(cardname): 
                annotations['card'][seqid] = cardname.split(',')
    
    return annotations



def parser_cazy(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:
    
            seqid, cazy_fam, cazy_genbank, cazy_species = map(str.strip, line.split('\t'))
            if clean(cazy_fam):
                annotations['cazy'][seqid] = cazy_fam.split(',')           
    
    return annotations


def parser_go(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:
    
            seqid, go_list = map(str.strip, line.split('\t'))
            if clean(go_list):
                annotations['go'][seqid] = go_list.split(',')
            
    return annotations


def parser_goslim(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:
            seqid, golist, godesclist, gocounter = map(str.strip, line.split('\t'))
            if clean(golist):
                annotations['goslim'][seqid] = golist.split(',')
    
    return annotations



def parser_pdb(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:
            seqid, pdbname = map(str.strip, line.split('\t'))
            if clean(pdbname):
                annotations['pdb'][seqid] = pdbname.split(',')
    
    return annotations


def parser_pfam(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:

            seqid, pfam_string = map(str.strip, line.rstrip('\n').split('\t'))

            pfam_array = []
            for pfam in pfam_string.split(','):
                #pfam_name, start, end = pfam.split('_')
                pfaminfo = pfam.split('_')
                pfam_name = '_'.join(pfaminfo[:-2])
                start = pfaminfo[-2]
                end = pfaminfo[-1]
                start = int(start)
                end = int(end)
                pfam_array.append([pfam_name, start, end])

            annotations['pfam'][seqid] = tuple(pfam_array)
            # Creates a signature of domain architecture based on domain names present in each seq
            annotations['pfam_arc'][seqid] = list([a[0] for a in pfam_array])
        
        return annotations



def parser_smart(file_path, annotations):

    with open(file_path, "r") as file:
        for line in file:

            seqid, smart_string = map(str.strip, line.rstrip('\n').split('\t'))

            smart_array = []
            for smart in smart_string.split(','):
                smart_name, start, end = smart.split('|')
                smart_name = smart_name.strip()
                start = int(start)
                end = int(end)
                smart_array.append([smart_name, start, end])

            annotations['smart'][seqid] = tuple(smart_array)
            # Creates a signature of domain architecture based on domain names present in each seq
            annotations['smart_arc'][seqid] = list([a[0] for a in smart_array])
        
        return annotations




           



