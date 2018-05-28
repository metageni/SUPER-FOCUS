# SUPER-FOCUS version 0.27 - This program is used to download and format the databases for the tool
import os,sys

#returns the path for a given program name
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def checkAligners():
    return [info for info in [which(x) for x in ["prerapsearch","diamond","makeblastdb"]] if info !=None]

def donwloadDBnFormat(aligners):
    print "Downloading DB"
    os.system("wget edwards.sdsu.edu/superfocus/downloads/db.zip")
    os.system("mv db.zip db")
    print "\nUncompressing DB"
    os.system("unzip db/db.zip") #uncompress db
    os.system("mv clusters/ db/")#mv db
    os.system("rm db/db.zip")#delete orignal downloaded file
    print "\nJoining files"
    [os.system("cat db/clusters/"+cluster+"/* > db/clusters/"+cluster+".fasta") for cluster in os.listdir("db/clusters/")]

    #Format database
    print "\nFormatting DB"
    for dbname in ["100","98","95","90"]:
        for aligner in aligners:
            if aligner == "prerapsearch":
                print "RAPSearch2: DB_"+dbname
                os.system("prerapsearch -d db/clusters/"+dbname+"_clusters.fasta -n db/static/rapsearch2/"+dbname+".db")
            elif aligner == "makeblastdb":
                print "blast: DB_"+dbname
                os.system("makeblastdb -in db/clusters/"+dbname+"_clusters.fasta -dbtype prot -out db/static/blast/"+dbname+".db -title "+dbname+".db")
            else:
                print "DIAMOND: DB_"+dbname
                os.system("diamond makedb --in  db/clusters/"+dbname+"_clusters.fasta --db db/static/diamond/"+dbname+".db")
    os.system("rm db/clusters/*.fasta")
try:
    aligners=[]
    temp_aligners=[x.split("/")[-1] for x in checkAligners() if x!="None"]
    
    given=[aligner.lower() for aligner in sys.argv[1:]]#which aligner the user wants to formart the database to
    c=0
    for i in sys.argv[1:]:
        if i.lower() not in ["rapsearch","blast","diamond"]:
            print i," is not a valid option of aligner; rapsearch, blast, or diamond are valid choices"
            c=1
    if c==0:
        if len(temp_aligners)>0:
            if "rapsearch" in given and "prerapsearch" in temp_aligners:aligners.append("prerapsearch")
            if "diamond" in given and "diamond" in temp_aligners:aligners.append("diamond")
            if "blast" in given and "makeblastdb" in temp_aligners:aligners.append("makeblastdb")

            os.system("rm db/clusters/ -r 2> /dev/null")#delete folder if exists
            donwloadDBnFormat(aligners)
            print "Done! Now you can run superfocus.py"
        else:
            print "None of the required aligners are installed [rapsearch, blast, or diamond]"
except:
    print "USAGE: python superfocus__downloadDB.py aligner_choice\nExample: python superfocus__downloadDB.py rapsearch blast diamond"
