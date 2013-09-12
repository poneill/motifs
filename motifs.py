"""Turn Prodoric into a usable Python module."""
#Code could stand some cleaning, but is usable.

import os,csv,sys
sys.path.append("lib/utils")
print sys.path
print os.getcwd()
from utils import *
from collections import namedtuple
BindingSiteRecord = namedtuple("BindingSiteRecord",["genome", "tf", "operon",
                                                    "strand", "sequence",
                                                    "startpos", "endpos",
                                                    "rflank", "lflank",
                                                    "regulation"])

class BindingSite(str):
    # Basically subclass string for TFBS attributes
    def __new__(cls,genome,tf,operon,strand,sequence,startpos,endpos,rflank,lflank,regulation):
        # YANETUT:
        # http://stackoverflow.com/questions/820742/how-to-bestow-string-ness-on-my-class
        newobj = str.__new__(cls,sequence)
        newobj.genome = genome
        newobj.tf = tf
        newobj.operon= operon
        newobj.strand = strand
        newobj.startpos = startpos
        newobj.endpos = endpos
        newobj.lflank = lflank
        newobj.rflank = rflank
        newobj.regulation = regulation
        return newobj

class Organism(object):
    pass

csv_names = [f for f in os.listdir('/home/poneill/motifs') #make robust later
        if f.endswith(".csv") and not f == "dump.csv"]

for csv_name in csv_names:
    org_name = "_".join(csv_name.split(" ")[:2])
    file_name = os.path.join("/home/poneill/motifs",csv_name)
    if org_name.endswith("."):
        org_name = org_name[:-1]
    #omit header
    records = map(BindingSiteRecord._make, csv.reader(open(file_name,'rb')))[1:]
    binding_sites = [BindingSite(record.genome,record.tf,record.operon,
                                 record.strand,record.sequence,record.startpos,
                                 record.endpos,record.rflank,record.lflank,
                                 record.regulation)
                     for record in records]
    tf_names = set(bs.tf for bs in binding_sites)
    exec("%s = Organism()" % org_name)
    for tf_name in tf_names:
        tf_binding_sites = [bs for bs in binding_sites if bs.tf == tf_name]
        setattr(eval(org_name),tf_name,tf_binding_sites)

    setattr(eval(org_name),"tfs",list(tf_names))
        
        #exec_string = ("%s.%s = %s" %
        #(org_name,tf_name,tf_binding_sites))
        # exec(exec_string)

        #This is perhaps the most subtle bug I have ever encountered.
        #I retain it here for posterity.  Never use exec when there's
        #another option.
        
        
print "loaded motifs"
