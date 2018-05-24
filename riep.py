#!/usr/bin/python

# script to iterate over different orientations and distances of the protein with respect
# to the membrane and call gromacs for a short MD to calculate the interaction energy

import sys, os
import getopt
import shutil
from string import replace
from subprocess import Popen, check_call, PIPE, CalledProcessError

# JRA: added further options to avoid hard-coding file names, box size, number of atoms and required mindist
def usage():
    print("""usage: {0} -h
       {0} <pitch> <roll> <yaw> -d <mindist> -x <centre_x> -y <centre_y> -z <centre_z>
           -a <box_dim_x> -b <box_dim_y> -c <box_dim_z>
           -p <protein_gro_file> -m <membrane_gro_file>
           -n <num_prot_memb_atoms>
           [ -s <exe_suffix> ]""".format(progname))
    return

exitcode = { 'ok' : 0, 'warning' : 1, 'critical' : 2, 'unknown' : 3 }

progname = os.path.basename(__file__)

# JRA: added d, bx, by, bz, prot, memb, natoms
try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hs:d:x:y:z:a:b:c:p:m:n:", ["help", "suffix="])
except getopt.GetoptError as err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(exitcode['critical'])

# JRA: added pm_dist, box_size, prot_name, memb_name, num_atoms
centre = [None, None, None]
box_size = [None, None, None]
suffix = ""
prot_name = ""
memb_name = ""
pm_dist = None
num_atoms = None
try:
    for opt, arg in opts:
        print(opt, arg)
        if opt in ("-h", "--help"):
            usage()
            sys.exit(exitcode['ok'])
        elif opt in ("-s", "--suffix"):
            suffix = arg
        elif opt == "-x":
            centre[0] = float(arg)
        elif opt == "-y":
            centre[1] = float(arg)
        elif opt == "-z":
            centre[2] = float(arg)
        elif opt == "-a":
            box_size[0] = float(arg)
        elif opt == "-b":
            box_size[1] = float(arg)
        elif opt == "-c":
            box_size[2] = float(arg)
        elif opt == "-p":
            prot_name = arg
        elif opt == "-m":
            memb_name = arg
        elif opt == "-n":
            num_atoms = int(arg)
        elif opt == "-d":
            pm_dist = float(arg)
        else:
            print("Failed")
            assert False, "unhandled option"
except ValueError:
    sys.stderr.write("""Centre coordinates must be floating-point numbers (or convertible to such)
""")
    sys.exit(exitcode['critical'])

# JRA: not sure I've got these values right...
if len(args) != 3 or len(centre) != 3 or len(box_size) != 3:
    usage()
    sys.exit(exitcode['critical'])

# Increment of angle of rotation from starting position (degrees)
try:
    pitch = int(args[0]) # Increment of rotation about X
    roll = int(args[1])  # Increment of rotation about Y
    yaw = int(args[2])   # Increment of rotation about Z
except ValueError:
    sys.stderr.write("""Angles of rotation must be integers (or convertible to integers)
""")
    sys.exit(exitcode['critical'])

[pitch, roll, yaw] = [angle % 360 for angle in [pitch, roll, yaw]]
rotstring = "_{0:03d}_{1:03d}_{2:03d}".format(pitch, roll, yaw)

# Get the location of GROMACS
try:
    gromacsbindir = os.path.join(os.environ['EBROOTGROMACS'], 'bin')
except KeyError:
    sys.stderr.write("""The environment variable EBROOTGROMACS is not defined.
Please load a GROMACS module, then try again.
""")
    sys.exit(exitcode['critical'])

# Make sure the suffix is correct
editconf_exe = os.path.join(gromacsbindir, "editconf{0}".format(suffix))
if not os.path.isfile(editconf_exe):
    suggested_suffix = None
    for gromacsprogram in os.listdir(gromacsbindir):
        # Assume only one "editconf"
        if gromacsprogram.startswith("editconf"):
            suggested_suffix = replace(gromacsprogram, "editconf", "", 1)
            break
    if suggested_suffix is None:
        sys.stderr.write("Error: editconf is missing -- something is badly wrong!\n")
    else:
        sys.stderr.write("{0}: no such file -- perhaps suffix should be {1}?\n".format(editconf_exe, suggested_suffix))
    sys.exit(exitcode['critical'])

# output directory for putting all the files that we generate into
outdir = prot_name+'_rot-ene'
try: 
    os.mkdir(outdir)
except:
    pass

# open a log file to see what's happening
log_filename = rotstring + ".log"
log = open(log_filename,'w')

# JRA: changed file names and box size to use program inputs
# JRA: what is meant by "big box"?
# protein coordinates only (in big box)
prot = prot_name+'.gro'
# membrane coordinates only (in big box)
mb = memb_name+'.gro'
# box size consistent across protein and membrane coordinates
end = [str(coord)+' ' for coord in box_size]

data_to_keep = []
# new name for rotated protein coordinates according to current iteration
rotprot = os.path.join(outdir, replace(prot,'.gro',rotstring+'.gro'))
# rotate and re-centre protein
# note that all arguments to check_call must be strings, not ints or floats
rotate_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, 'editconf'), suffix) ] 
rotate_cmd.extend(["-f", prot])
rotate_cmd.extend(["-rotate", str(pitch), str(roll), str(yaw)])
rotate_cmd.extend(["-center"])
rotate_cmd.extend([str(coord) for coord in centre])
rotate_cmd.extend(["-o", rotprot ])
log.write("Running rotate command: {0}\n".format(rotate_cmd))
check_call(rotate_cmd)

# new name for rotated protein coordinates + membrane
rotprotmb = os.path.join(outdir, replace(prot,'.gro',rotstring+'mb.gro'))
# manually concatenate rotated protein coordinates with membrane: genbox doesn't work
# read in rotated protein coordinates
fi = open(rotprot,'r')
li = fi.readlines()
fi.close()
fm = open(mb,'r')
lm = fm.readlines()
fm.close()
# open output file for combined protein + membrane coordinates
fo = open(rotprotmb,'w')
# JRA: changed to use number of atoms from program inputs
# write header inclusive of atom count in protein + membrane
init = 'Rotated protein + membrane;\n'+str(num_atoms)+'\n'
fo.write(init)
# write rotated protein coordinates
for i in range(2,len(li)-1):
    fo.write(li[i])
# write membrane coordinates
for i in range(2,len(lm)-1):
    fo.write(lm[i])
# write box size
fo.writelines(end)
fo.close()

# output file for mindist
mdfile = os.path.join(outdir,"{0}_mindist{1}.xvg".format(prot_name,rotstring))
# calculate mininum distance between protein and membrane
# it asks for groups:
#Select a group: 1
#Selected 1: 'Protein'
#Select a group: 12
#Selected 12: 'Other'
# group selection hard-coded to 1 and 12:
while True:
    dist_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, 'g_mindist'), suffix) ] 
    dist_cmd.extend(['-f', rotprotmb])
    dist_cmd.extend(['-n', "{0}.ndx".format(prot_name)])
    dist_cmd.extend(['-od', mdfile])

    log.write("Running distance command: {0}\n".format(dist_cmd))
    proc = Popen(dist_cmd, stdin=PIPE)
    proc.communicate('1\n12\n')
    retcode = proc.wait()
    if retcode != 0:
        raise CalledProcessError(returncode=retcode, cmd=dist_cmd)

# extract mindist from file
    fd = open(mdfile,'r')
    # hard-coded to standard mindist file format
    mindist = float(fd.readlines()[19].split()[1])
    fd.close()
    log.write('Mindist = '+str(mindist)+'\n\n')

    # JRA: changed to use pm_dist from program inputs rather than have 0.50 hard-coded
    # move protein coordinates so that minimum distance between protein and membrane is exactly pm_dist nm
    # move protein away if it overlaps with membrane or is too close, and closer if it's too far away
    if mindist > pm_dist : transdist = round(mindist - pm_dist,3)
    # NOTE: may not round as expected so may not quite shift enough
    # JRA: What do you mean "may not round as expected"?
    else : transdist = -1*round(pm_dist - mindist,3)
    log.write('Transdist = '+str(transdist)+'\n\n')

    # translate to current working distance (only translate in z)
    transprot = "{0}_translated.gro".format(rotstring)
    translate_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, 'editconf'), suffix) ] 
    translate_cmd.extend(["-f", rotprot])
    translate_cmd.extend(["-translate", "0", "0", str(transdist)])
    translate_cmd.extend(["-o", transprot ])
    log.write("Running translate command: {0}\n".format(translate_cmd))
    check_call(translate_cmd)

    # replace rotated protein coordinates with translated protein coordinates
    os.rename(transprot, rotprot)
 
    # combine coordinates (overwriting previous)
    # read in translated protein coordinates
    fi = open(rotprot,'r')
    li = fi.readlines()
    fi.close()
    fm = open(mb,'r')
    lm = fm.readlines()
    fm.close()
    # open output file for combined protein + membrane coordinates (overwriting previous)
    fo = open(rotprotmb,'w')
    # JRA: changed to use number of atoms from program input
    # write header, inclusive of atom count in protein + membrane
    init = 'Rotated protein + membrane;\n'+str(num_atoms)+'\n'
    fo.write(init)
    # write protein coordinates
    for i in range(2,len(li)-1):
        fo.write(li[i])
    # write membrane coordinates
    for i in range(2,len(lm)-1):
        fo.write(lm[i])

    # write box size
    fo.writelines(end)
    fo.close()

    # re-calculate protein-membrane distance
    log.write("Running distance command: {0}\n".format(dist_cmd))
    proc = Popen(dist_cmd, stdin=PIPE)
    proc.communicate('1\n12\n')
    retcode = proc.wait()
    if retcode != 0:
        raise CalledProcessError(returncode=retcode, cmd=dist_cmd)
    #check_call(dist_cmd)

    # extract new mindist from file
    fd = open(mdfile,'r')
    mindist = float(fd.readlines()[19].split()[1])
    fd.close()
    log.write('New mindist = '+str(mindist)+'\n\n')

    if -0.01 < transdist < 0.01 : break

# make separate copy of protein + membrane topology file for each rotation
shutil.copy("{0}.top".format(prot_name), "{0}{1}.top".format(prot_name,rotstring))

# Solvate each rotated protein + membrane system
solvate_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, "genbox"), suffix) ]
solvate_cmd.extend([ '-cp', rotprotmb ])
solvate_cmd.extend([ '-cs', 'spc216.gro' ])
solvate_cmd.extend([ '-p', "{0}{1}.top".format(prot_name,rotstring) ])
solvate_cmd.extend([ '-o', "{0}_solv{1}.gro".format(prot_name,rotstring) ])
log.write("Running solvate command: {0}\n".format(solvate_cmd))
check_call(solvate_cmd)

# Energy minimise solvated system
min_grompp_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, 'grompp'), suffix) ]
min_grompp_cmd.extend([ '-f', "minim.mdp" ])
min_grompp_cmd.extend([ '-c', "{0}_solv{1}.gro".format(prot_name,rotstring) ])
min_grompp_cmd.extend([ '-p', "{0}{1}.top".format(prot_name,rotstring) ])
min_grompp_cmd.extend([ '-o', "{0}_em{1}.tpr".format(prot_name,rotstring) ])

log.write("Preprocessing energy minimisation: {0}\n".format(min_grompp_cmd))
check_call(min_grompp_cmd)

# Run energy minimisation
min_mdrun_cmd = [ "srun" ]
min_mdrun_cmd.extend([ "{0}{1}".format(os.path.join(gromacsbindir, 'mdrun'), suffix) ])
min_mdrun_cmd.extend([ '-deffnm', "{0}_em{1}".format(prot_name,rotstring) ])
log.write("Running energy minimisation: {0}\n".format(min_mdrun_cmd))
check_call(min_mdrun_cmd)
    
# prepare for 25-step MD
md_grompp_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, 'grompp'), suffix) ]
md_grompp_cmd.extend([ '-f', "final.mdp" ])
md_grompp_cmd.extend([ '-c', "{0}_em{1}.gro".format(prot_name,rotstring) ])
md_grompp_cmd.extend([ '-p', "{0}{1}.top".format(prot_name,rotstring) ])
md_grompp_cmd.extend([ '-o', "{0}{1}_energy".format(prot_name,rotstring) ])
log.write("Preprocessing MD simulation: {0}".format(md_grompp_cmd))
check_call(md_grompp_cmd)

# run 25-step MD
md_mdrun_cmd = [ "srun" ]
md_mdrun_cmd.extend([ "{0}{1}".format(os.path.join(gromacsbindir, 'mdrun'), suffix) ])
md_mdrun_cmd.extend([ '-deffnm', "{0}{1}_energy".format(prot_name,rotstring) ])
log.write("Running MD simulation: {0}".format(md_mdrun_cmd))
check_call(md_mdrun_cmd)

# calculate interaction energy (vdw as well as electrostatic)
# hard-coded for standard order of energies
# 47  Coul-SR:Protein-Other
# 48  LJ-SR:Protein-Other
# 49  LJ-LR:Protein-Other
md_energy_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, "g_energy"), suffix) ]
md_energy_cmd.extend([ '-f', "{0}{1}_energy.edr".format(prot_name, rotstring) ])
md_energy_cmd.extend([ '-s', "{0}{1}_energy.tpr".format(prot_name, rotstring) ])
md_energy_cmd.extend([ '-o', "{0}{1}_energy.xvg".format(prot_name, rotstring) ])
log.write("Computing energies: {0}".format(md_energy_cmd))
proc = Popen(md_energy_cmd, stdin=PIPE)
proc.communicate('47\n48\n49\n')
retcode = proc.wait()
if retcode != 0:
    raise CalledProcessError(returncode=retcode, cmd=md_energy_cmd)

# if desired, select residues within cutoff distance (1.4 nm)
#select_cmd = [ "{0}{1}".format(os.path.join(gromacsbindir, "g_select"), suffix) ]
#select_cmd.extend([ '-f', "{0}{1}_energy.gro".format(prot_name, rotstring) ])
#select_cmd.extend([ '-s', "{0}{1}_energy.tpr".format(prot_name, rotstring) ])
#select_cmd.extend([ '-select', "'\"reslect\" residues within 1.4 of group other'" ])
#select_cmd.extend([ '-oi', "{0}{1}.dat".format(prot_name, rotstring) ])
#log.write("Selecting residues: {0}".format(select_cmd))
#check_call(select_cmd)

with open("{0}{1}_energy.xvg".format(prot_name, rotstring),'r') as fh:
    important_line = fh.readlines()[-1].split()
data_to_keep.append('%d,%d,%d\t%s\t%s\t%s'%(pitch,roll,yaw,important_line[1],important_line[2],important_line[3]))

with open('outputfile_'+str(rotstring)+'.txt','w') as fh:
    for i in data_to_keep: fh.write(i+'\n')

# delete unnecessary files
to_delete = [f for f in os.listdir(".") if "{0}_energy.trr".format(rotstring) in f]
for file in to_delete:
    os.remove(file)

to_delete = [f for f in os.listdir(".") if "_em{0}".format(rotstring) in f]
for file in to_delete:
    os.remove(file)

to_delete = [f for f in os.listdir(".") if "_solv{0}".format(rotstring) in f]
for file in to_delete:
    os.remove(file)

to_delete = [f for f in os.listdir(".") if "mdout" in f]
for file in to_delete:
    os.remove(file)

sys.exit(exitcode['ok'])
