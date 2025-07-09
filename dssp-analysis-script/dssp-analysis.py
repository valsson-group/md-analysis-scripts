#! /usr/bin/env python 

import mdtraj as md
import numpy as np
import argparse 
import time

parser = argparse.ArgumentParser(conflict_handler="resolve")

parser.add_argument("--pdb-file",             type=str, required=True)
parser.add_argument("--trajectory-files",     type=str, required=True,  nargs="*")
parser.add_argument("--residue-set",          type=str, required=False, nargs="*")
parser.add_argument("--output-file",          type=str, required=False, default="dssp.data")
parser.add_argument("--stride",               type=int, required=False, default=1)
parser.add_argument("--output-chains",        type=str, required=False, nargs="*")
parser.add_argument("--simplifed",            action="store_true", required=False)

args = parser.parse_args()

traj_files=args.trajectory_files
pdb_file=args.pdb_file 
simplified=args.simplifed 
dssp_fn = args.output_file
stride=args.stride
# dssp_fracton_fn = "dssp.fraction.data"

LoopsAndIrregular_Code="L"
if not simplified:
    AssignmentCodes="#! SET dssp_type full\n"
else:
    AssignmentCodes="#! SET dssp_type simplifed\n"

AssignmentCodes+="""#
# DSSP assignment codes 
# \"H\" : Alpha helix
# \"B\" : Residue in isolated beta-bridge
# \"E\" : Extended strand, participates in beta ladder
# \"G\" : 3-helix (3/10 helix)
# \"I\" : 5 helix (pi helix)
# \"T\" : hydrogen bonded turn
# \"S\" : bend
# \"{:1s}\" : Loops and irregular elements
# """.format(LoopsAndIrregular_Code)

if simplified:
    AssignmentCodes+="""
# Using simplifed DSSP assignment codes
# \"H\" : Helix. Either of the \"H\", \"G\", or \"I\" codes.
# \"E\" : Strand. Either of the \"E\", or \"B\" codes.
# \"C\" : Coil. Either of the \"T\", \"S\" or \"{:1s}\" codes
# """.format(LoopsAndIrregular_Code)


chain_labels = ["a","b","c","d","e","f","g"]

f_dssp = []
# To make sure that the `total_sasa_fn` file
# is started clean for the run
if args.output_chains:
    for c in args.output_chains:
        dssp_fn_chain = ".chain-{:s}.".format(c.lower()).join(dssp_fn.rsplit(".",1 ))
        f_tmp = open(dssp_fn_chain,"w")
        f_tmp.close()
        f_dssp.append(open(dssp_fn_chain,"a"))
else:
    f_tmp = open(dssp_fn,"w")
    f_tmp.close()
    f_dssp.append(open(dssp_fn,"a"))

if args.residue_set:
    dssp_fn_residue_set = ".residue-set.".join(dssp_fn.rsplit(".",1 ))
    f_dssp_residue_set = open(dssp_fn_residue_set,"w")
    f_dssp_residue_set.close()
    f_dssp_residue_set = open(dssp_fn_residue_set,"a")

number_of_traj_files=len(traj_files)
counter=0
for trj in traj_files:
    start_time_loop = time.perf_counter() 
    counter+=1
    print(" {0:02d}/{1:02d}: {2:s}".format(counter,number_of_traj_files,trj), flush=True)

    traj = md.load_xtc(trj, 
                       pdb_file,
                       stride=stride)
    dssp = md.compute_dssp(traj,
                           simplified=simplified)
    if not simplified:
        dssp[dssp == " "] = LoopsAndIrregular_Code
    time_traj = traj.time
    time_traj.shape=(-1,1)

    header_time = "#! FIELDS time "
    residues_str = [] 
    for i,ch in  enumerate(traj.topology.chains):
        ch = chain_labels[i]
        for r in traj.topology.chain(i).residues:
            res_str = "{0:1s}-{1:s}".format(ch,str(r).lower())
            residues_str.append(res_str)
    #

    residue_set_id = []
    residue_set_str = []
    if args.residue_set:
        for res in args.residue_set:
            try:
                index = residues_str.index(res.lower())
                residue_set_id.append(index)
                residue_set_str.append(res.lower())
            except ValueError as ve:
                print("Warning: cannot find {:s} among the residue labels, ignored".format(res))

    res_indices = []
    res_strings = []
    if args.output_chains:
        for ch in args.output_chains:
            id_tmp = [i for i in range(len(residues_str)) if residues_str[i][0:1] == ch.lower()]
            str_tmp = [r for r in residues_str if r[0:1] == ch.lower()]
            if len(id_tmp)==0: 
                print("Error: cannot find chain {:s} among the chains".format(ch))
                exit()
            res_indices.append(id_tmp)
            res_strings.append(str_tmp)
    else:
        id_tmp = [i for i in range(len(residues_str))]
        res_indices.append(id_tmp)
        res_strings.append(residues_str)

    for f,ids,res in zip(f_dssp,res_indices,res_strings):
        np.savetxt(f,
                   X = np.concatenate((time_traj,dssp[:,ids]),axis=1),
                   header=header_time + " ".join(res) + "\n"+AssignmentCodes,
                   comments="",
                   fmt="%s")
        f.flush()

    if args.residue_set:
        np.savetxt(f_dssp_residue_set,
                X = np.concatenate((time_traj,dssp[:,residue_set_id]),axis=1),
                   header=header_time + " ".join(residue_set_str) + "\n"+AssignmentCodes,
                   comments="",
                   fmt="%s")
        f_dssp_residue_set.flush()
    end_time_loop = time.perf_counter() 
    print(" - Elapsed time: {:0.1f} seconds".format(end_time_loop-start_time_loop), flush=True)

for f in f_dssp:
    f.close()
if args.residue_set: f_dssp_residue_set.close()
