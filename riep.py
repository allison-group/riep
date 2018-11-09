#!/usr/bin/python

# Iterates over different orientations and distances of a protein with respect to a membrane
# Calls GROMACS to run a short MD simulation to calculate the interaction energy

import sys, os
import getopt
import shutil
from string import replace
from subprocess import Popen, check_call, PIPE, CalledProcessError

def usage():
    print("""usage: {0} -h
       {0} <taitBryan_x> <taitBryan_y> <taitBryan_z> -d <mindist> 
           -x <box_dim_x> -y <box_dim_y> -z <box_dim_z>
           -p <protein_gro_file> -m <membrane_gro_file>
           -n <num_prot_memb_atoms>
           -r <exe_prefix>
           [ -s <exe_suffix> ]""".format(progname))
    return

exitcode = { 'ok' : 0, 'warning' : 1, 'critical' : 2, 'unknown' : 3 }

progname = os.path.basename(__file__)

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hs:d:x:y:z:p:m:n:r:", ["help", "suffix=", "prefix="])
except getopt.GetoptError as err:
    # print help information and exit:
    print str(err)
    usage()
    sys.exit(exitcode['critical'])

pm_dist = None
box_size = [None, None, None]
prot_name = ""
memb_name = ""
num_atoms = None
prefix = ""
suffix = ""
try:
    for opt, arg in opts:
        print(opt, arg)
        if opt in ("-h", "--help"):
            usage()
            sys.exit(exitcode['ok'])
        elif opt == "-d":
            pm_dist = float(arg)
        elif opt == "-x":
            box_size[0] = float(arg)
        elif opt == "-y":
            box_size[1] = float(arg)
        elif opt == "-z":
            box_size[2] = float(arg)
        elif opt == "-p":
            prot_name = arg
        elif opt == "-m":
            memb_name = arg
        elif opt == "-n":
            num_atoms = int(arg)
        elif opt in ("-r", "--prefix"):
            prefix = arg
        elif opt in ("-s", "--suffix"):
            suffix = arg
        else:
            print("Failed")
            assert False, "unhandled option"
except ValueError:
    sys.stderr.write("""Protein-membrane distance and centre and box coordinates must be floating-point numbers (or convertible to such)
""")
    sys.exit(exitcode['critical'])

if len(args) != 3 or len(box_size) != 3:
    usage()
    sys.exit(exitcode['critical'])

# Set rotation angle (wrt starting position; degrees)
try:
    taitBryan_x = int(args[0]) # Increment of rotation about X
    taitBryan_y = int(args[1])  # Increment of rotation about Y
    taitBryan_z = int(args[2])   # Increment of rotation about Z
except ValueError:
    sys.stderr.write("""Angles of rotation must be integers (or convertible to integers)
""")
    sys.exit(exitcode['critical'])

[taitBryan_x, taitBryan_y, taitBryan_z] = [angle % 360 for angle in [taitBryan_x, taitBryan_y, taitBryan_z]]
rotstring = "_{0:03d}_{1:03d}_{2:03d}".format(taitBryan_x, taitBryan_y, taitBryan_z)

# Set box centre
centre = [box_size[0]/2., box_size[1]/2., box_size[2]/2.]

# Get the location of GROMACS
try:
    gromacsbindir = os.path.join(os.environ['EBROOTGROMACS'], 'bin')
except KeyError:
    sys.stderr.write("""The environment variable EBROOTGROMACS is not defined.
Please load a GROMACS module, then try again.
""")
    sys.exit(exitcode['critical'])

# Set up program names
mindist_exe = ""
editconf_exe = ""
make_ndx_exe = ""
solvate_exe = ""
grompp_exe = ""
mdrun_exe = ""
energy_exe = ""
select_exe = ""
if (prefix == 'g_'):
    mindist_exe = os.path.join(gromacsbindir, "{0}mindist{1}".format(prefix,suffix))
    editconf_exe = os.path.join(gromacsbindir, "editconf{0}".format(suffix))
    make_ndx_exe = os.path.join(gromacsbindir, "make_ndx{0}".format(suffix))
    solvate_exe = os.path.join(gromacsbindir, "genbox{0}".format(suffix))
    grompp_exe = os.path.join(gromacsbindir, "grompp{0}".format(suffix))
    mdrun_exe = os.path.join(gromacsbindir, "mdrun{0}".format(suffix))
    energy_exe = os.path.join(gromacsbindir, "{0}energy{1}".format(prefix,suffix))
    select_exe = os.path.join(gromacsbindir, "{0}select{1}".format(prefix,suffix))
elif (prefix == 'gmx'):
    mindist_exe = os.path.join(gromacsbindir, "{0}{1} mindist".format(prefix,suffix))
    editconf_exe = os.path.join(gromacsbindir, "{0}{1} editconf".format(prefix,suffix))
    make_ndx_exe = os.path.join(gromacsbindir, "{0}{1} make_ndx".format(prefix,suffix))
    solvate_exe = os.path.join(gromacsbindir, "{0}{1} solvate".format(prefix,suffix))
    grompp_exe = os.path.join(gromacsbindir, "{0}{1} grompp".format(prefix,suffix))
    mdrun_exe = os.path.join(gromacsbindir, "{0}{1} mdrun".format(prefix,suffix))
    energy_exe = os.path.join(gromacsbindir, "{0}{1} energy".format(prefix,suffix))
    select_exe = os.path.join(gromacsbindir, "{0}{1} select".format(prefix,suffix))
else:
    sys.stderr.write("""Prefix must be supplied.""")
    sys.exit(exitcode['critical'])

# Make sure the prefix and suffix are correct
# Testing suffix with editconf
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

# Testing prefix with mindist
if not os.path.isfile(mindist_exe):
    mindist_exists = False
    for gromacsprogram in os.listdir(gromacsbindir):
        # Assume only one "mindist"
        if (gromacsprogram.find("mindist") != -1):
            mindist_exists = True
            break
    if mindist_exists is False:
        sys.stderr.write("Error: mindist is missing -- something is badly wrong!\n")
    else:
        sys.stderr.write("{0}: no such file -- perhaps prefix is wrong?\n".format(mindist_exe))
    sys.exit(exitcode['critical'])


# Make output directory for putting all the files that we generate into
outdir = prot_name+'_rot-ene'
try: 
    os.mkdir(outdir)
except:
    pass

# Open a log file to see what's happening
log_filename = rotstring + ".log"
log = open(log_filename,'w')

# protein coordinates only (box must be large enough for entire system)
prot = prot_name+'.gro'
# membrane coordinates only (box must be large enough for entire system)
mb = memb_name+'.gro'
# box size consistent across protein and membrane coordinates
end = [str(coord)+' ' for coord in box_size]

data_to_keep = []
# new name for rotated protein coordinates according to current iteration
rotprot = os.path.join(outdir, replace(prot,'.gro',rotstring+'.gro'))
# Rotate and re-centre protein
# note that all arguments to check_call must be strings, not ints or floats
rotate_cmd = [ editconf_exe ] 
rotate_cmd.extend(["-f", prot])
rotate_cmd.extend(["-rotate", str(taitBryan_x), str(taitBryan_y), str(taitBryan_z)])
rotate_cmd.extend(["-center"])
rotate_cmd.extend([str(coord) for coord in centre])
rotate_cmd.extend(["-o", rotprot ])
log.write("Running rotate command: {0}\n".format(rotate_cmd))
check_call(rotate_cmd)

# new name for rotated protein coordinates + membrane
rotprotmb = os.path.join(outdir, replace(prot,'.gro',rotstring+'_'+memb_name+'.gro'))
# Manually concatenate rotated protein coordinates with membrane: genbox doesn't work
# Read in rotated protein coordinates
fi = open(rotprot,'r')
li = fi.readlines()
fi.close()
fm = open(mb,'r')
lm = fm.readlines()
fm.close()
# Open output file for combined protein + membrane coordinates
fo = open(rotprotmb,'w')
# Write header inclusive of atom count in protein + membrane
init = 'Rotated protein + membrane;\n '+str(num_atoms)+'\n'
fo.write(init)
# Write rotated protein coordinates
for i in range(2,len(li)-1):
    fo.write(li[i])
# Write membrane coordinates
for i in range(2,len(lm)-1):
    fo.write(lm[i])
# Write box size
fo.writelines(end)
fo.close()

# Create index file
ndx_name = prot_name+'_'+memb_name
if not os.path.isfile("{0}.ndx".format(ndx_name)):
    ndx_cmd = [ make_ndx_exe ]
    ndx_cmd.extend(['-f', rotprotmb])
    ndx_cmd.extend(['-o', "{0}.ndx".format(ndx_name)])

    log.write("Making index file: {0}\n".format(ndx_cmd))
    proc = Popen(ndx_cmd, stdin=PIPE)
    proc.communicate('q\n')
    retcode = proc.wait()
    if retcode != 0:
        raise CalledProcessError(returncode=retcode, cmd=ndx_cmd)

# output file for mindist
mdfile = os.path.join(outdir,"{0}_mindist{1}.xvg".format(prot_name,rotstring))
# Calculate mininum distance between protein and membrane
# NOTE: group selection hard-coded to 1 ('Protein') and 12 ('Other'):
while True:
    dist_cmd = [ mindist_exe ] 
    dist_cmd.extend(['-f', rotprotmb])
    dist_cmd.extend(['-n', "{0}.ndx".format(ndx_name)])
    dist_cmd.extend(['-od', mdfile])

    log.write("Running distance command: {0}\n".format(dist_cmd))
    proc = Popen(dist_cmd, stdin=PIPE)
    proc.communicate('1\n12\n')
    retcode = proc.wait()
    if retcode != 0:
        raise CalledProcessError(returncode=retcode, cmd=dist_cmd)

    # Extract mindist from file
    fd = open(mdfile,'r')
    # NOTE: hard-coded to standard mindist file format
    mindist = float(fd.readlines()[19].split()[1])
    fd.close()
    log.write('Mindist = '+str(mindist)+'\n\n')

    # Move protein coordinates so that minimum distance between protein and membrane is exactly pm_dist nm
    # Move protein away if it overlaps with membrane or is too close, and closer if it's too far away
    if mindist > pm_dist : transdist = round(mindist - pm_dist,3)
    # NOTE: due to rounding, may not quite shift enough but iterations should fix this
    else : transdist = -1*round(pm_dist - mindist,3)
    log.write('Transdist = '+str(transdist)+'\n\n')

    # Translate to current working distance (only translate in z)
    transprot = "{0}_translated.gro".format(rotstring)
    translate_cmd = [ editconf_exe ] 
    translate_cmd.extend(["-f", rotprot])
    translate_cmd.extend(["-translate", "0", "0", str(transdist)])
    translate_cmd.extend(["-o", transprot ])
    log.write("Running translate command: {0}\n".format(translate_cmd))
    check_call(translate_cmd)

    # Replace rotated protein coordinates with translated protein coordinates
    os.rename(transprot, rotprot)
 
    # Combine coordinates (overwriting previous)
    # Read in translated protein coordinates
    fi = open(rotprot,'r')
    li = fi.readlines()
    fi.close()
    fm = open(mb,'r')
    lm = fm.readlines()
    fm.close()
    # Open output file for combined protein + membrane coordinates (overwriting previous)
    fo = open(rotprotmb,'w')
    # Write header, inclusive of atom count in protein + membrane
    init = 'Rotated protein + membrane;\n'+str(num_atoms)+'\n'
    fo.write(init)
    # Write protein coordinates
    for i in range(2,len(li)-1):
        fo.write(li[i])
    # Write membrane coordinates
    for i in range(2,len(lm)-1):
        fo.write(lm[i])

    # Write box size
    fo.writelines(end)
    fo.close()

    # Re-calculate protein-membrane distance
    log.write("Running distance command: {0}\n".format(dist_cmd))
    proc = Popen(dist_cmd, stdin=PIPE)
    proc.communicate('1\n12\n')
    retcode = proc.wait()
    if retcode != 0:
        raise CalledProcessError(returncode=retcode, cmd=dist_cmd)
    #check_call(dist_cmd)

    # Extract new mindist from file
    fd = open(mdfile,'r')
    mindist = float(fd.readlines()[19].split()[1])
    fd.close()
    log.write('New mindist = '+str(mindist)+'\n\n')

    if -0.01 < transdist < 0.01 : break

# Make separate copy of protein + membrane topology file for this rotation
shutil.copy("{0}.top".format(prot_name), "{0}{1}.top".format(prot_name,rotstring))

# Solvate this rotated protein + membrane system
solvate_cmd = [ solvate_exe ]
solvate_cmd.extend([ '-cp', rotprotmb ])
solvate_cmd.extend([ '-cs', 'spc216.gro' ])
solvate_cmd.extend([ '-p', "{0}{1}.top".format(prot_name,rotstring) ])
solvate_cmd.extend([ '-o', "{0}_solv{1}.gro".format(prot_name,rotstring) ])
log.write("Running solvate command: {0}\n".format(solvate_cmd))
check_call(solvate_cmd)

# Energy minimise solvated system
min_grompp_cmd = [ grompp_exe ]
min_grompp_cmd.extend([ '-f', "minim.mdp" ])
min_grompp_cmd.extend([ '-c', "{0}_solv{1}.gro".format(prot_name,rotstring) ])
min_grompp_cmd.extend([ '-p', "{0}{1}.top".format(prot_name,rotstring) ])
min_grompp_cmd.extend([ '-o', "{0}_em{1}.tpr".format(prot_name,rotstring) ])

log.write("Preprocessing energy minimisation: {0}\n".format(min_grompp_cmd))
check_call(min_grompp_cmd)

# Run energy minimisation
# NOTE: uncomment next two lines and comment out following line to run on cluster with slurm
#min_mdrun_cmd = [ "srun" ]
#min_mdrun_cmd.extend([ mdrun_exe ])
min_mdrun_cmd = [ mdrun_exe ]
min_mdrun_cmd.extend([ '-deffnm', "{0}_em{1}".format(prot_name,rotstring) ])
log.write("Running energy minimisation: {0}\n".format(min_mdrun_cmd))
check_call(min_mdrun_cmd)
    
# Prepare for 25-step MD
md_grompp_cmd = [ grompp_exe ]
md_grompp_cmd.extend([ '-f', "final.mdp" ])
md_grompp_cmd.extend([ '-c', "{0}_em{1}.gro".format(prot_name,rotstring) ])
md_grompp_cmd.extend([ '-p', "{0}{1}.top".format(prot_name,rotstring) ])
md_grompp_cmd.extend([ '-o', "{0}{1}_energy".format(prot_name,rotstring) ])
log.write("Preprocessing MD simulation: {0}".format(md_grompp_cmd))
check_call(md_grompp_cmd)

# Run 25-step MD
# NOTE: uncomment next two lines and comment out following line to run on cluster with slurm
#md_mdrun_cmd = [ "srun" ]
#md_mdrun_cmd.extend([ mdrun_exe ])
md_mdrun_cmd = [ mdrun_exe ]
md_mdrun_cmd.extend([ '-deffnm', "{0}{1}_energy".format(prot_name,rotstring) ])
log.write("Running MD simulation: {0}".format(md_mdrun_cmd))
check_call(md_mdrun_cmd)

# Calculate interaction energy (vdw and electrostatic)
# NOTE: hard-coded for standard order of energies
# 47  Coul-SR:Protein-Other
# 48  LJ-SR:Protein-Other
# 49  LJ-LR:Protein-Other
md_energy_cmd = [ energy_exe ]
md_energy_cmd.extend([ '-f', "{0}{1}_energy.edr".format(prot_name, rotstring) ])
md_energy_cmd.extend([ '-s', "{0}{1}_energy.tpr".format(prot_name, rotstring) ])
md_energy_cmd.extend([ '-o', "{0}{1}_energy.xvg".format(prot_name, rotstring) ])
log.write("Computing energies: {0}".format(md_energy_cmd))
proc = Popen(md_energy_cmd, stdin=PIPE)
proc.communicate('47\n48\n49\n')
retcode = proc.wait()
if retcode != 0:
    raise CalledProcessError(returncode=retcode, cmd=md_energy_cmd)

# OPTIONAL: select residues within cutoff distance (1.4 nm)
# Ideally only for one or a few energetically optimal rotations
#select_cmd = [ select_exe ]
#select_cmd.extend([ '-f', "{0}{1}_energy.gro".format(prot_name, rotstring) ])
#select_cmd.extend([ '-s', "{0}{1}_energy.tpr".format(prot_name, rotstring) ])
#select_cmd.extend([ '-select', "'\"reslect\" residues within 1.4 of group other'" ])
#select_cmd.extend([ '-oi', "{0}{1}.dat".format(prot_name, rotstring) ])
#log.write("Selecting residues: {0}".format(select_cmd))
#check_call(select_cmd)

with open("{0}{1}_energy.xvg".format(prot_name, rotstring),'r') as fh:
    important_line = fh.readlines()[-1].split()
data_to_keep.append('%d,%d,%d\t%s\t%s\t%s'%(taitBryan_x,taitBryan_y,taitBryan_z,important_line[1],important_line[2],important_line[3]))

with open('outputfile_'+str(rotstring)+'.txt','w') as fh:
    for i in data_to_keep: fh.write(i+'\n')

# Delete unnecessary files
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
