import sys,os,subprocess
from subprocess import Popen, PIPE, STDOUT
import pexpect;
from Bio import SeqIO;
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord;

def is_right_column_to_keep(column):
    """
    Determines whether a column should be kept or deleted
    """
    total_unique=0;
    if( column.upper().count("-")>=1 ):

        return False;
    for i in range(len(column)):
        if column.upper().count( column[i] )==1:
            total_unique=total_unique+1;
            #print column[i];

    #print total_unique, "/", len(column);
    return (total_unique/len(column))<0.5;



def fasta_to_phylip(infilename,outfilename):
    """
    convert the original alignment to noise reduced phylip alignment
    """
    handle = open(infilename, "rU");
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"));
    handle.close()


    i=0;
    new_records=[];

    column_dict={};
    column="";
    first_iteration=True; # at first Iteration ,initialise column dictionary
    id_list=[];
    for k in record_dict.keys():

        sequence=record_dict[k].seq;
        id_list.append(k);  # save the sequence ID's for secure ordering
        if first_iteration: # at first Iteration ,initialise column list
            #print "At first Iteration ...";
            first_iteration=False;
            for j in range(len(sequence)):
                column_dict[j]="";

        #print "ID: " ,record_dict[k].id , "Len: ", len(sequence);
        #print "Dict Length: ", len( column_dict.keys() );
        for j in  range(len(column_dict.keys()) ):
            column_dict[j]=column_dict[j]+sequence[j];

        i=i+1;


    seq_dict={};
    all_column_deleted=True; # Flag to test whether all column deleted or not
    first_iteration=True; # at first Iteration ,initialise sequence dict
    for j in  range(len(column_dict.keys()) ): # iterate over all the column
        current_column=column_dict[j];

        if first_iteration: # at first Iteration ,initialise sequence dict
            first_iteration=False;
            for id in id_list:
                seq_dict[id]="";


        if(is_right_column_to_keep(current_column) ): # is it a good column ?
            all_column_deleted=False;

            for  m  in  range(len(id_list)): # insert one character to each of sequence from current column
                seq_dict[id_list[m]]=seq_dict[ id_list[m] ]+ current_column[m] ;



    if all_column_deleted:
        print "Error: All columns are deleted from  ",infilename;
        sys.exit();
    for id in id_list:
        if( len(seq_dict[id])==0 ):
            print "Error: ",infilename;
            return ;
        record=SeqRecord( Seq(seq_dict[id]),id ,'','')
        new_records.append(record);


    output_handle = open(outfilename, "w");
    SeqIO.write(new_records, output_handle, "phylip");
    output_handle.close();

    return ;

def phylip_protdist(infilename,replica):
    """
    Calculates the distance matrix from for a given phylip formated alignment file
    """

    prog = 'phylip protdist'

    responses = 'rs.txt'

    FH = open(responses,'w')

    current_dir = os.getcwd()

    FH.write (current_dir + '/' + infilename + '\n')

    if pexpect.os.path.exists(current_dir + '/infile'):
        #print 'infile exists';
        pexpect.os.remove(current_dir+'/infile');

    if pexpect.os.path.exists(current_dir + '/outfile'):
        #print 'outfile exists';
        FH.write('R\n')

    if (replica>1):
        FH.write('M\n')
        FH.write('D\n')
        FH.write(str(replica)+'\n')
    FH.write('Y\n')
    FH.close()

    cmd = prog
    cmd += ' < ' + responses
    cmd += ' >  screenout '


    #print cmd ;

    p = subprocess.Popen(cmd, shell=True)

    pid,ecode = os.waitpid(p.pid, 0)

    outfilename="newout.dist"

    pexpect.os.rename('outfile', outfilename);

    if pexpect.os.path.exists(current_dir+'/outfile'):
        #print 'removing outfile';
        pexpect.os.remove(current_dir+'/outfile');

    if pexpect.os.path.exists(current_dir+'/'+infilename):
        #print 'removing outfile';
        pexpect.os.remove(current_dir+'/'+infilename);

    return outfilename;

def phylip_neighbor(infilename,script_dir,result_group,treename,replica):
    """
    Construct phylogenic tree from distance matrix
    """

    prog ='phylip neighbor'

    responses = 'rs.txt'

    FH = open(responses,'w')

    current_dir = os.getcwd()

    #print current_dir;

    FH.write (current_dir + '/' + infilename + '\n')

    if pexpect.os.path.exists(current_dir + '/infile'):
        #print 'infile exists';
        pexpect.os.remove(current_dir+'/infile');

    if pexpect.os.path.exists(current_dir+'/outfile'):
        #print 'outfile exists';
        FH.write('R\n')

    if (replica>1):
        FH.write('M\n')

        FH.write(str(replica)+'\n')
        FH.write('4339\n')


    FH.write('Y\n')


    if pexpect.os.path.exists(current_dir + '/outtree'):
        #print 'outfile exists';
        FH.write('R\n')


    FH.close()

    cmd = prog
    cmd += ' < ' + responses
    cmd += ' >  screenout '




    p = subprocess.Popen(cmd, shell=True)

    pid,ecode = os.waitpid(p.pid, 0)

    if pexpect.os.path.exists(current_dir+'/outfile'):
        pexpect.os.remove(current_dir+'/outfile');

    if pexpect.os.path.exists(current_dir+'/rs.txt'):
        pexpect.os.remove(current_dir+'/rs.txt');


    if pexpect.os.path.exists(current_dir+'/outtree'):

        pexpect.os.rename(current_dir+'/outtree', treename+'.tree');

    if pexpect.os.path.exists(current_dir+'/screenout'):
        pexpect.os.remove(current_dir+'/screenout');
    if pexpect.os.path.exists(current_dir+'/'+infilename):
        pexpect.os.remove(current_dir+'/'+infilename);


    if pexpect.os.path.exists(current_dir+'/newtree'):
        pexpect.os.remove(current_dir+'/newtree');

    return "tree";
