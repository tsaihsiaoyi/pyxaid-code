# ***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
# ***********************************************************/

import os
import math

# Module contains:
# def calculate(energy_prefix,energy_suffix,dip_prefix,dip_suffix,isnap,fsnap,opt,scl1,scl2,outfile,HOMO,minE,maxE,dE)
# def ham_map(prefix,isnap,fsnap,suffix,opt,scl,outfile)
# def ham_map1(energy_prefix,energy_suffix,prefix,isnap,fsnap,suffix,opt,scl1,scl2,outfile)


def calculate(energy_prefix, energy_suffix, dip_prefix, dip_suffix, isnap, fsnap, opt, scl1, scl2, outfile, HOMO, minE, maxE, dE):
    # Compute average value of the Hamiltonian matrix elements in energy axis:
    # energy_prefix - prefix of the files containing the energy of the states (diagonal terms)
    # energy_suffix - --//--
    # dip_prefix - is the prefix of the series of filenames containing transition dipole moment
    # dip_suffix - is the suffix of the series of filenames containing transition dipole moment
    # isnap - index of initial file
    # fsnap - index of final file
    # opt - option for averaging
    # scl1 - scaling factor for energy scale
    # scl2 - scaling factor for output quantity scale
    # outfile - is the name of the file, where the map will be written
    # HOMO - is the index (starting from 1) of the HOMO orbital
    # minE, maxE, dE determine the range and the resolution of the absorption spectrum

    filename = energy_prefix + str(isnap) + energy_suffix
    f = open(filename, "r")
    a = f.readline()
    f.close()
    sz = len(a.split())
    print("Map size is %d by %d" % (sz, sz))

    HOMO = HOMO - 1  # now it has the meaning of the index rather than id

    # ========= Prepare storage ============
    M = []
    i = 0
    while i < sz:
        m = []
        j = 0
        while j < sz:
            m.append(0.0)
            j = j + 1
        M.append(m)
        i = i + 1

    exE = []
    exI = []
    e = minE
    npt = 0
    while e < (maxE+dE):
        exE.append(e)
        exI.append(0.0)
        e = e + dE
        npt = npt + 1

    # ========== Read in data ==============
    count = 0.0

    for k in range(isnap, fsnap+1):
        filename = dip_prefix + str(k) + dip_suffix
        e_filename = energy_prefix + str(k) + energy_suffix
        print("Reading files ", filename, e_filename)

        if(os.path.exists(filename) and os.path.exists(e_filename)):

            # Read energies of the states
            f1 = open(e_filename, "r")
            A1 = f1.readlines()
            f1.close()

            E = []
            i = 0
            while i < sz:
                tmp1 = A1[i].split()
                y = scl1*float(tmp1[i])
                E.append(y)
                i = i + 1

            # Read transition dipole moments
            f = open(filename, "r")
            A = f.readlines()
            f.close()
            i = 0
            while i < sz:
                tmp = A[i].split()
                j = 0
                while j < sz:
                    x = float(tmp[j])
                    if opt == 0:
                        M[i][j] = x
                    elif opt == 1:
                        M[i][j] = math.fabs(x)
                    elif opt == 2:
                        M[i][j] = x*x

                    j = j + 1
                i = i + 1

            # Now compute transition intensities

            i = HOMO
            while i > 0:
                j = HOMO+1
                while j < sz:
                    # ====== Now we are at the point to consider i->j transition
                    de = E[j] - E[i]
                    if(de > minE and de < maxE):
                        indx = int((de - minE)/dE)
                        exI[indx] = exI[indx] + M[i][j]
                        count = count + 1.0

                    j = j + 1
                i = i-1

    # ========== Average and print out results ===============
    out = open(outfile, "w")

    exI[0] = 0.0   # Do not print dipole moments
    i = 0
    while i < npt:
        if(opt == 0 or opt == 1):
            exI[i] = scl2*exI[i] / count
        elif opt == 2:
            exI[i] = scl2 * math.sqrt(exI[i]) / count

        out.write("%f  %f \n" % (exE[i], exI[i]))
#        line = str(exE[i]) + "  "+str(exI[i]) + "\n"
#        out.write(line)
        i = i + 1

    out.close()

    return [exE, exI]


# Example of usage
#opt = 2
# scl1 = 13.60569253 # Ry to eV
# scl2 = 0.12 # scaling for transition dipole moments to give eV*fs units
#HOMO = 12
#minE = 0.0
#maxE = 4.0
#dE = 0.1

# calculate("res/0_Ham_","_re","res/0_Hprime_","x_re",0,6,opt,scl1,scl2,"ab_spectrx.dat",HOMO,minE,maxE,dE)


def ham_map(prefix, isnap, fsnap, suffix, opt, scl, outfile):
    # Compute average value of the Hamiltonian matrix elements in energy axis:
    # prefix - is the prefix of the series of filenames containing required matrix elements
    # isnap - index of initial file
    # fsnap - index of final file
    # suffix - similarly to prefix
    # opt - option for averaging
    # scl - scaling factor. e.g. to convert units
    # outfile - is the name of the file, where the map will be written

    filename = prefix + str(isnap) + suffix
    f = open(filename, "r")
    a = f.readline()
    f.close()
    sz = len(a.split())
    print("Map size is %d by %d" % (sz, sz))

    # ========= Prepare storage ============
    M = []
    i = 0
    while i < sz:
        m = []
        j = 0
        while j < sz:
            m.append(0.0)
            j = j + 1
        M.append(m)
        i = i + 1

    # ========== Read in data ==============
    count = 0.0

    for k in range(isnap, fsnap+1):
        filename = prefix + str(k) + suffix
        print("Reading file ", filename)

        if(os.path.exists(filename)):

            f = open(filename, "r")
            A = f.readlines()
            f.close()

            i = 0
            while i < sz:
                tmp = A[i].split()
                j = 0
                while j < sz:
                    x = float(tmp[j])
                    if opt == 0:
                        M[i][j] = M[i][j] + x
                    elif opt == 1:
                        M[i][j] = M[i][j] + math.fabs(x)
                    elif opt == 2:
                        M[i][j] = M[i][j] + x*x

                    j = j + 1
                i = i + 1

            count = count + 1.0

    # ========== Average and print results ===============
    out = open(outfile, "w")

    i = 0
    while i < sz:
        j = 0
        while j < sz:
            # ====== Average ======
            if(opt == 0 or opt == 1):
                M[i][j] = scl * M[i][j] / count
            elif opt == 2:
                M[i][j] = scl * math.sqrt(M[i][j] / count)

            # ======= Print =======
            if(i == j):
                # because diagonal and off-diagonal elements are of different orders of magnitude
                M[i][j] = 0.0
            out.write("%d  %d  %f\n" % (i, j, M[i][j]))
            j = j + 1
        out.write("\n")
        i = i + 1

    return 1


def ham_map1(energy_prefix, energy_suffix, prefix, isnap, fsnap, suffix, opt, scl1, scl2, outfile):
    # Compute average value of the Hamiltonian matrix elements in energy axis:
    # energy_prefix - prefix of the files containing the energy of the states (diagonal terms)
    # energy_suffix - --//--
    # prefix - is the prefix of the series of filenames containing required matrix elements
    # isnap - index of initial file
    # fsnap - index of final file
    # suffix - similarly to prefix
    # opt - option for averaging
    # scl1 - scaling factor for energy scale
    # scl2 - scaling factor for output quantity scale
    # outfile - is the name of the file, where the map will be written

    filename = prefix + str(isnap) + suffix
    f = open(filename, "r")
    a = f.readline()
    f.close()
    sz = len(a.split())
    print("Map size is %d by %d" % (sz, sz))

    # ========= Prepare storage ============
    M = []
    E = []
    i = 0
    while i < sz:
        m = []
        j = 0
        while j < sz:
            m.append(0.0)
            j = j + 1
        M.append(m)
        E.append(0.0)
        i = i + 1

    # ========== Read in data ==============
    count = 0.0

    for k in range(isnap, fsnap+1):
        filename = prefix + str(k) + suffix
        e_filename = energy_prefix + str(k) + energy_suffix
        print("Reading files ", filename, e_filename)

        if(os.path.exists(filename) and os.path.exists(e_filename)):

            f = open(filename, "r")
            A = f.readlines()
            f.close()

            f1 = open(e_filename, "r")
            A1 = f1.readlines()
            f1.close()

            i = 0
            while i < sz:
                tmp = A[i].split()
                tmp1 = A1[i].split()
                y = float(tmp1[i])
                j = 0
                while j < sz:
                    x = float(tmp[j])
                    if opt == 0:
                        M[i][j] = M[i][j] + x
                    elif opt == 1:
                        M[i][j] = M[i][j] + math.fabs(x)
                    elif opt == 2:
                        M[i][j] = M[i][j] + x*x

                    j = j + 1
                E[i] = E[i] + y
#                print(i,y)
                i = i + 1

            count = count + 1.0

    # ========== Average and print results ===============
    out = open(outfile, "w")

    i = 0
    while i < sz:
        E[i] = scl1 * E[i] / count
        i = i + 1

    i = 0
    while i < sz:
        j = 0
        while j < sz:
            # ====== Average ======
            if(opt == 0 or opt == 1):
                M[i][j] = scl2 * M[i][j] / count
            elif opt == 2:
                M[i][j] = scl2 * math.sqrt(M[i][j] / count)

            # ======= Print =======
            if(i == j):
                # because diagonal and off-diagonal elements are of different orders of magnitude
                M[i][j] = 0.0
            out.write("%f  %f  %f\n" % (E[i], E[j], M[i][j]))
            j = j + 1
        out.write("\n")
        i = i + 1

    out.close()

    return 1


# Example of usage
#opt = 2
# scl1 = 13.60569253 # Ry to eV
# scl2 = 0.12 # scaling for transition dipole moments to give eV*fs units
#
#ham_map("res/0_Ham_",   0,6,"_re" ,opt,scl1,"ave_Ham_re.dat")
#ham_map("res/0_Ham_",   0,6,"_im" ,opt,scl1,"ave_Ham_im.dat")
# ham_map("res/0_Hprime_",0,6,"x_re",opt,scl2,"ave_Hprime_x_re.dat")
# ham_map("res/0_Hprime_",0,6,"y_re",opt,scl2,"ave_Hprime_y_re.dat")
# ham_map("res/0_Hprime_",0,6,"z_re",opt,scl2,"ave_Hprime_z_re.dat")
#
#ham_map1("res/0_Ham_","_re","res/0_Ham_",   0,6,"_re" ,opt,scl1,scl1,"1ave_Ham_re.dat")
#ham_map1("res/0_Ham_","_re","res/0_Ham_",   0,6,"_im" ,opt,scl1,scl1,"1ave_Ham_im.dat")
# ham_map1("res/0_Ham_","_re","res/0_Hprime_",0,6,"x_re",opt,scl1,scl2,"1ave_Hprime_x_re.dat")
# ham_map1("res/0_Ham_","_re","res/0_Hprime_",0,6,"y_re",opt,scl1,scl2,"1ave_Hprime_y_re.dat")
# ham_map1("res/0_Ham_","_re","res/0_Hprime_",0,6,"z_re",opt,scl1,scl2,"1ave_Hprime_z_re.dat")
#
