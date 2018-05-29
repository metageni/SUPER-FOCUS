#!/usr/bin/env python2
"""SUPER-FOCUS version 0.29. This program is used to download and format the databases for SUPER-FOCUS."""
import os, sys

# returns the path for a given program name
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
    return [
        info
        for info in [which(x) for x in ["prerapsearch", "diamond", "makeblastdb"]]
        if info != None
    ]


def donwloadDBnFormat(aligners):
    print "Downloading DB"
    os.system("wget edwards.sdsu.edu/superfocus/downloads/db.zip")
    os.system("mv db.zip db")
    print "\nUncompressing DB"
    os.system("unzip db/db.zip")  # uncompress db
    os.system("mv clusters/ db/")  # mv db
    os.system("rm db/db.zip")  # delete orignal downloaded file
    print "\nJoining files"
    [
        os.system(
            "cat db/clusters/" + cluster + "/* > db/clusters/" + cluster + ".fasta"
        )
        for cluster in os.listdir("db/clusters/")
    ]

    # Format database
    print "\nFormatting DB"
    cluster_identities = ["100", "98", "95", "90"]
    for dbname in cluster_identities:
        for aligner in aligners:
            if aligner == "prerapsearch":
                print "RAPSearch2: DB_" + dbname
                os.system(
                    "prerapsearch -d db/clusters/"
                    + dbname
                    + "_clusters.fasta -n db/static/rapsearch2/"
                    + dbname
                    + ".db"
                )
            elif aligner == "makeblastdb":
                print "blast: DB_" + dbname
                os.system(
                    "makeblastdb -in db/clusters/"
                    + dbname
                    + "_clusters.fasta -dbtype prot -out db/static/blast/"
                    + dbname
                    + ".db -title "
                    + dbname
                    + ".db"
                )
            else:
                print "DIAMOND: DB_" + dbname
                os.system(
                    "diamond makedb --in  db/clusters/"
                    + dbname
                    + "_clusters.fasta --db db/static/diamond/"
                    + dbname
                    + ".db"
                )
    os.system("rm db/clusters/*.fasta")


def main():
    valid_aligners = {"rapsearch", "blast", "diamond"}
    aligner_db_creators = {os.path.basename(x) for x in checkAligners() if x != "None"}
    if not aligner_db_creators:
        print "None of the required aligners are installed", list(valid_aligners)
        sys.exit(1)

    # Parse which aligner(s) the user wants to format the database to
    requested_aligners = {aligner.lower() for aligner in sys.argv[1:]}
    for aligner in requested_aligners:
        if aligner not in allowed_aligners:
            print i, " is not a valid option of aligner; 'rapsearch', 'blast', or 'diamond' are valid choices"
            sys.exit(1)
    if "rapsearch" in given and "prerapsearch" in aligner_db_creators:
        aligners.append("prerapsearch")
    if "diamond" in given and "diamond" in aligner_db_creators:
        aligners.append("diamond")
    if "blast" in given and "makeblastdb" in aligner_db_creators:
        aligners.append("makeblastdb")

    os.system("rm db/clusters/ -r 2> /dev/null")  # delete folder if exists
    donwloadDBnFormat(aligners)
    print "Done! Now you can run superfocus.py"
    sys.exit()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print __doc__
        print "USAGE: python superfocus__downloadDB.py aligner_choice"
        print "Example: python superfocus__downloadDB.py rapsearch blast diamond"
        sys.exit(1)
    main()
