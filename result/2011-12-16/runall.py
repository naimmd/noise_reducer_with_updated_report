#!/usr/bin/env python
import sys,os,subprocess
from subprocess import Popen, PIPE, STDOUT
import pexpect;
from Bio import SeqIO;
import tempfile;
from Bio.SeqRecord import SeqRecord
import BootStrapLib;
import dendropy;
from dendropy import treecalc;



replica =1 ;
script_dir = os.getcwd();

data_dir = "../../data/2011-12-16/";
result_dir="../../result/2011-12-16/";


os.chdir(data_dir); # Go to data directory first

command="ls"
datagroups=pexpect.os.popen(command,'r').read();

print "Frequency of reference tree recovery";
print "------------------------------------"
print "Data Group    ","        Before    ","    After    " ; #

# lists to keep track of average difference
avg_distance_before=[];
avg_distance_after=[];

#Indicator whether we really have some data group or not
no_data_group=True;
data_group=[];

ref_org=[]; # Reference to Original tree difference
ref_red=[]; #Reference to Reduced tree difference

for eachgroup in datagroups.split():
    #empty the list for each data group
    del ref_org[:];
    del ref_red[:];
    no_data_group=False; #That means you have some data group to work with;

    os.chdir(eachgroup); # go inside each data group of data directory
    data_group_directory=os.getcwd();
    command="ls *.msl"
    allalignments=pexpect.os.popen(command,'r').read();

    outfilename="temp.phy";


    if os.path.exists(script_dir+"/"+eachgroup)==False:
        os.mkdir(script_dir+"/"+eachgroup);
    else: # Romove al existing files
        allfiles=pexpect.os.popen("ls  "+script_dir+"/"+eachgroup ).read();
        for each_file in  allfiles.split():
            pexpect.os.remove(script_dir+"/"+eachgroup+"/"+each_file);



    #taxa = dendropy.TaxonSet();

    ref_tree_name=pexpect.os.popen("ls *.tree",'r').read().split()[0];

    try:
        ref_tree = dendropy.Tree.get_from_path(ref_tree_name, 'newick'); #,taxon_set=taxa
    except Exception:
        print "\n  Something wrong with Reference tree or the file is empty";
        sys.exit();

    is_first_alignment=True;
    first_record_count=0;
    record_count=0;
    for each_alignment in allalignments.split():


        #Write Original and noise reduced alignments to file in phylip format


        try:
            record_count=SeqIO.convert(each_alignment,"fasta",script_dir+"/"+eachgroup+"/"+outfilename, "phylip"); #Original alignment

            if is_first_alignment==True:
                first_record_count=record_count;
                is_first_alignment=False;

        except ValueError:
            print "\n *Error in Input file, Finishing execution* \n TODO: Sequences have different length or input File is empty or not properly formated for sequences :",eachgroup,"/", each_alignment;
            print " \n";
            sys.exit();

        if is_first_alignment==False:
            if first_record_count !=record_count:
                print "\n*Error Detected , Stopping  execution* \n TODO: We Expect same number of records in all alignment file \n";
                sys.exit();

        if record_count<=1:
            print "\n*Error Detected , Stopping  execution* \n Are you Serious !! only one or no  record in your alignment file named : ",eachgroup,"/", each_alignment;
            sys.exit();


        BootStrapLib.fasta_to_phylip(each_alignment, script_dir+"/"+eachgroup+"/"+"R"+outfilename); # Noise reduced alignment

        os.chdir(script_dir+"/"+eachgroup); # Go to result group directory

        #TODO; Get original and noise reduced tree

        #First get distance matrix
        #Second ; From distance matrix to  tree


        filename=BootStrapLib.phylip_protdist(outfilename,replica); #Protdist for Original alignment
        tree=BootStrapLib.phylip_neighbor(filename,script_dir,eachgroup,each_alignment,replica); #neighbor for Original alignment

        try:
            original_tree=dendropy.Tree.get_from_path(each_alignment+".tree", 'newick'); #,taxon_set=taxa  tree  Original alignment
        except IndexError:
            print "\n*Error Detected , Stopping  execution* \n TODO: Tree from original sequence got to be mysterious for dendropy \n";
            sys.exit();





        filename=BootStrapLib.phylip_protdist("R"+outfilename,replica); #Protdist for Original alignment
        tree=BootStrapLib.phylip_neighbor(filename,script_dir,eachgroup,"R"+each_alignment,replica); #neighbor for Original alignment

        noise_reduced_tree=dendropy.Tree.get_from_path("R"+each_alignment+".tree", 'newick'); #,taxon_set=taxa  tree  Original alignment



        # Save the difference from reference tree to original and noise reduced tree
        ref_org.append( ref_tree.symmetric_difference(original_tree) );
        ref_red.append( ref_tree.symmetric_difference(noise_reduced_tree));



        os.chdir(data_group_directory); # get back to data group directory


    #show how often they are recovered
    if eachgroup[0]=='s':
        print "",
    print eachgroup , "       ",ref_org.count(0),"            ",ref_red.count(0);
    data_group.append(eachgroup);
    avg_distance_before.append( float(sum(ref_org))/len(ref_org)) ;
    avg_distance_after.append(float(sum(ref_red))/len(ref_red) );



    os.chdir(".."); # back to parent data directory

print "------------------------------------"
print "Average distance ";
print "------------------"
print "Data Group    ","        Before    ","    After    " ; #
for j in range(len(data_group)):
    if data_group[j][0]=='s':
        print "",
    print data_group[j],"       ",avg_distance_before[j] ,"            ",avg_distance_after[j];

#print "ref_org:",ref_org;
#print "ref_red:",ref_red;

if(no_data_group==True):
    print "\n Well ; what are we suppose to do without any data group!!\n";