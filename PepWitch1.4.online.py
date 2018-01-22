print("Loading Modules...")
import numpy as np
import glob
import os.path
import time
import datetime
import itertools
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib_venn import venn2
from math import log10
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests
from statistics import mean
import bisect
from sys import argv
import shutil
import sys

def main():
    np.set_printoptions(precision=2)
    Disregard = int(argv[2])
    MinSpc = float(argv[1])
    SpectralFraction = float(argv[5])
    SameSamePValue = float(argv[3])
    WorkingDirectory = os.path.dirname(os.path.realpath(__file__))
    TTest = False
    UniqueA = False
    UniqueB = False
    Engine = argv[7]
    CSVInputName = argv[8:]
    Control = CSVInputName[0]
    if len(CSVInputName) > 1:
        Wild = CSVInputName[1:]
    DesiredReplicateNumber = int(argv[4])

    def GetInputName():
        filenames = 0
        control = 0
        wild = 0
        replicates = []

        print("Welcome to PyWitch!")
        print("You can either import a single state for analysis or perform a two-state comparison. Right now, only GPM table *.csv files are supported.")
        print("For a two-state comparison, please enter the control state first, then the wild type.")
        print("######################################")

        while filenames == 0:
            filenames = input("Provide the name of the file stems that you want to import, separated by a space: ").split()

            for name in filenames:
                if not glob.glob(str(WorkingDirectory) + "/" + str(name) + "*.csv"):
                    print("######################################")
                    print("Error #1: " + str(name) + " not found. The program will look for every filename provided a common stem.")
                    print("So for example, if you type in 'gr', the program will pull in both 'grape.csv' and 'green.csv', etc")
                    filenames = 0

            if not filenames:
                filenames = 0

            try:
                if len(filenames) != len(set(filenames)):
                    print("######################################")
                    print("Error #5: Files found, but you've entered the same name more than once. Why. That's silly. Don't do that.")
                    filenames = 0
            except TypeError:
                continue

            try:
                if len(filenames) >= 3:
                    multistatetypelock = 0
                    while multistatetypelock == 0:
                        print("Do you want to conduct a {control vs multistate} analysis or an {anova-style} analysis")
                        try:
                            multistatetype = int(input("[1] for multistate, [2] for anova, [3] to return to the beginning: "))
                            if int(multistatetype) == 3:
                                filenames = 0
                            if int(multistatetype)  == 2:
                                filenames = 0
                            if int(multistatetype)  == 1:
                                print("Which state is your control: ")
                                for i in range(len(filenames)):
                                    print(str(filenames[i]) + " = " + str(i+1))
                                controlindex = 0
                                while controlindex == 0:
                                    controlstate = input()
                                    try:
                                        if 1 <= int(controlstate) <= len(filenames):
                                            control = filenames[int(controlstate)-1]
                                            wild = [filenames[x] for x in range(len(filenames)) if x != int(controlstate)-1]
                                            print(control)
                                            print(wild)
                                            multistatetypelock = 1
                                            controlindex = 1
                                        else:
                                            print("######################################")
                                            print("Error #8: Please choose an appropriate control from the list.")
                                            controlindex = 0
                                    except TypeError:
                                        print("Something is going wrong.")
                                        controlindex = 0
                            elif int(multistatetype) != 1 or 2 or 3:
                                print("######################################")
                                print("Error #9: Invalid number entered. Try again.")
                        except ValueError:
                            print("######################################")
                            print("Error #6: Invalid input. Please try again")
                            multistatetypelock = 0
                        if filenames == 0:
                            multistatetypelock = 1
            except TypeError:
                print("error")
                filenames = 0

        for x in range(len(filenames)):
            replicates.append(len(glob.glob(str(WorkingDirectory) + "/" + str(filenames[x]) + "*.csv")))

        if replicates[:] != replicates[:]:
            print("######################################")
            print("Error #2: You have specified two sample states with different variable numbers: " + ''.join(str(replicates)).strip('[]'))
            print("The program will now exit. Please ensure that your replicate numbers are the same before retrying.")
            exit()

        replicatecarry = replicates[0]

        if len(filenames) == 2:
            control = filenames[0]
            wild = filenames[1]

        return filenames, control, wild, replicatecarry

    # CSVInputName, Control, Wild, DesiredReplicateNumber = GetInputName()

    DayandTime = os.path.dirname(os.path.realpath(__file__)) + "/output/" + str(datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S')) + '-' + '-'.join(CSVInputName) + "-" + str(MinSpc)
    os.makedirs(DayandTime)

    starttime = time.perf_counter()

    def ImportProteinFiles(inputnames): ## This class imports the files into our workflow and makes a list of the replicates.
        a = inputnames
        TotalDataIDsL = []
        UniqueNames = []
        TotalDataIDsD = []

        for name in a:

            if Engine == "GPM":
                try:
                    total = [np.genfromtxt(x, delimiter = ',', skip_header = 1, usecols = [0,2,5,3],  dtype = ('U5', 'float32','float32','uint32')) for x in glob.glob(str(argv[6]) + "/" + str(name) + "*.csv")] # Imports Protein Name, Spectral Count, Mol. Weight, score
                except ValueError:
                    print("######################################")
                    print("Error #4: Your CSV file is empty or incompatible. Currently, only GPM table *.csv outputs are supported." )
                    exit()
            elif Engine == "PD":
                try:
                    total = [np.genfromtxt(x, delimiter = ',', skip_header = 1, usecols = [0,7,9,2],  dtype = ('U5', 'float32','float32','uint32')) for x in glob.glob(str(argv[6]) + "/" + str(name) + "*.csv")] # Imports Protein Name, PSM, Mol. Weight, score
                except ValueError:
                    print("######################################")
                    print("Error #4: Your CSV file is empty or incompatible. Currently, only GPM table *.csv outputs are supported." )
                    exit()
            replicatenumber = len(total)
            # print("This is the number of replicates I found for " + str(name) + ": " + str(replicatenumber))

            flatlist = [] # List of the total imports list but without numpy. Numpy and strings don't mix particularly well.
            for i in range(replicatenumber):
                flatlist.append(total[i].tolist()) # Collapse numpy list so we can make sense of all the data; shape and replicate format still intact.

            uniquenamedict = {}
            totaldict = [ dict((row[0], row[1:]) for row in t) for t in flatlist ]

            for data in totaldict:
                uniquenamedict.update(data)

            TotalDataIDsL.append(flatlist)
            UniqueNames.append(uniquenamedict)
            TotalDataIDsD.append(totaldict)

        return TotalDataIDsL, UniqueNames, TotalDataIDsD

    TotalDataIDsList, UniqueDictionary, TotalDataDictionary = ImportProteinFiles(CSVInputName)

    def CollectUniqueIdsAndSpectra(udict, tdict, spec):   # This will use the unique protein list and match all IDs from the total list. NSAF, SpC, AA, MW.

        print("Finding Unique IDs...")

        spectralcount = []
        filteredvalues = []

        for y in range(len(udict)):

            spectralcountdict = {}
            for name in udict[y]:
                spectralcountdict[name] = [x[0] for x in (data.get(name, [spec]) for data in tdict[y])], [b for b in (next(x for x in (data.get(name, [spec])[1:2] for data in tdict[y]) if x))], [i[2] for i in (data.get(name, [spec, 0, 1]) for data in tdict[y])]
            spectralcount.append(spectralcountdict)

            filteredvaluesdict = {} # This is Scrappy, basically. Conducts as minimum spectral count and disregard replicate drop function
            for name in udict[y]:
                if ([x[0] for x in (data.get(name, [0]) for data in tdict[y])].count(0) <= Disregard) and (sum(x[0] for x in (data.get(name, [0]) for data in tdict[y])) >= MinSpc):
                    filteredvaluesdict[name] = [x[0] for x in (data.get(name, [spec]) for data in tdict[y])], [b for b in (next(x for x in (data.get(name, [spec])[1:2] for data in tdict[y]) if x))], [i[2] for i in (data.get(name, [spec, 0, 1]) for data in tdict[y])]
            filteredvalues.append(filteredvaluesdict)

        return spectralcount, filteredvalues

    UniqueIDsAndSpCDictionary, FilteredIDsAndSpCDictionary = CollectUniqueIdsAndSpectra(UniqueDictionary, TotalDataDictionary, SpectralFraction)

    def GetNSAF(a, vanilla, filtered):

        print("Calculating NSAFs...")

        logNSAFDictList = []
        logNSAFDictListFiltered = []

        for y in range(len(a)): # For every test state [control, wild]

            NSAF = []
            SAF = []
            nloge = []
            SAFBuffer = []
            NSAFBuffer = []
            nlogebuffer = []
            SumSPC = []
            CumSaf = []
            logNSAFDict = {}


            SumMol = (sum(float((1000*k[1][0])/110) for v, k in a[y].items()))

            for q in range(len(vanilla[y])): # For every replicate in state y
                SumSPC.append(sum(k[0][q] for v,k in a[y].items() if k[0][q] != SpectralFraction))

            for value in a[y]:
                SAFBuffer2 = []
                for q in range(len(vanilla[y])):  # For every replicate in state y
                    SAFBuffer2.append( ((a[y][value][0][q]) / float((1000*a[y][value][1][0])/110)) / (SumSPC[q] / SumMol) )
                SAFBuffer.append([SAFBuffer2, value])
            SAF.append(SAFBuffer)

            for q in range(len(vanilla[y])): # For every replicate in state y
                CumSaf.append(sum((o[0][q] for o in SAF[0])))  # <-- Calculating the total SAF values for making NSAF

            for t in range(len(vanilla[y])):
                nlogebuffer.append(sum(k[2][t] for v, k in a[y].items()))

            for p in a[y]:
                nlogebuffer2 = []
                for i in range(len(vanilla[y])):
                    nlogevalue = a[y][p][2][i] / nlogebuffer[i]
                    nlogebuffer2.append(nlogevalue)
                nloge.append(nlogebuffer2)

            for r in range(len(SAF[0])):
                NSAFBuffer2 = []
                for q in range(len(vanilla[y])):  # For every replicate in state y
                    NSAFBuffer2.append( log10(((SAF[0][r][0][q]) / CumSaf[q]))) # <-- Calculates the logNSAF value
                NSAFBuffer.append([SAF[0][r][1], NSAFBuffer2])
            NSAF.append(NSAFBuffer)

            for index,item in enumerate(a[y]):
                logNSAFDict[item] = [ p for p in a[y][item][0] ] , a[y][item][1] , [ u for u in NSAF[0][index][1] ], [ p for p in a[y][item][2]] , [e for e in nloge[index] ]
            logNSAFDictList.append(logNSAFDict)

        for y in range(len(filtered)): # For every test state [control, wild]

            NSAF = []
            SAF = []
            nloge = []
            SAFBuffer = []
            NSAFBuffer = []
            nlogebuffer = []
            SumSPC = []
            CumSaf = []
            logNSAFDict = {}

            SumMol = (sum(float((1000*k[1][0])/110) for v, k in filtered[y].items()))

            for q in range(len(vanilla[y])): # For every replicate in state y
                SumSPC.append(sum(k[0][q] for v,k in filtered[y].items()))

            for value in filtered[y]:
                SAFBuffer2 = []
                for q in range(len(vanilla[y])):  # For every replicate in state y
                    SAFBuffer2.append( ((filtered[y][value][0][q]) / float((1000*filtered[y][value][1][0])/110)) / (SumSPC[q] / SumMol) )
                SAFBuffer.append([SAFBuffer2, value])
            SAF.append(SAFBuffer)

            for q in range(len(vanilla[y])): # For every replicate in state y
                CumSaf.append(sum((o[0][q] for o in SAF[0])))  # <-- Calculating the total SAF values for making NSAF

            for r in range(len(SAF[0])): # For every SAF value
                NSAFBuffer2 = []
                for q in range(len(vanilla[y])):  # For every replicate in state y
                    NSAFBuffer2.append( log10(((SAF[0][r][0][q]) / CumSaf[q]))) # <-- Calculates the logNSAF value
                NSAFBuffer.append([SAF[0][r][1], NSAFBuffer2])
            NSAF.append(NSAFBuffer)

            for t in range(len(vanilla[y])):
                nlogebuffer.append(sum(k[2][t] for v, k in filtered[y].items()))

            for p in filtered[y]:
                nlogebuffer2 = []
                for i in range(len(vanilla[y])):
                    nlogevalue = filtered[y][p][2][i] / nlogebuffer[i]
                    nlogebuffer2.append(nlogevalue)
                nloge.append(nlogebuffer2)

            for index,item in enumerate(filtered[y]):
                logNSAFDict[item] = [ p for p in filtered[y][item][0] ] , filtered[y][item][1] , [ u for u in NSAF[0][index][1] ], [ p for p in filtered[y][item][2] ], [e for e in nloge[index] ]
            logNSAFDictListFiltered.append(logNSAFDict)

        return(logNSAFDictList, logNSAFDictListFiltered)

    UniqueIDsSpcNSAFDict, UniqueIDsSpcNSAFFilteredDict = GetNSAF(UniqueIDsAndSpCDictionary, TotalDataIDsList, FilteredIDsAndSpCDictionary)

    def DoTTest(names, filter, control, wild, Q):

        if Q == -1:
            print("Now conducting TTests for all proteins...")
        elif Q != -1:
            print("Now conducting Adjusted QTests for all proteins...")

        if len(filter) > 2:

            uniqueida = []
            uniqueidb = []
            ttestlist = []
            qtestlist = []

            for state in wild:

                qtestlista = []
                ttestdict = {}
                qtestdict = {}

                c = names.index(control)
                w = names.index(state)

                setstateControl = set(filter[c])
                setstateWild = set(filter[w])

                if Q == -1:
                    for name in setstateControl.intersection(setstateWild):
                        FC = mean(filter[c][name][2]) / mean(filter[w][name][2])
                        ttestdict[name] = stats.ttest_ind(filter[c][name][2], filter[w][name][2])[1], filter[c][name], filter[w][name], FC
                    ttestlist.append(ttestdict)
                elif Q != -1:
                    for name in setstateControl.intersection(setstateWild):
                        ttestdict[name] = stats.ttest_ind(filter[c][name][2], filter[w][name][2])[1],
                    tvals = [i[0] for i in ttestdict.values()]
                    qtestlista.append(multipletests(tvals, method='fdr_bh')[1])
                    g = 0
                    for name in setstateControl.intersection(setstateWild):
                        FC = mean(filter[c][name][2]) / mean(filter[w][name][2])
                        qtestdict[name] = qtestlista[0][g], filter[c][name], filter[w][name], FC
                        g += 1
                    qtestlist.append(qtestdict)

                uniqueA = {}
                uniqueB = {}

                for name in (setstateControl | setstateWild):
                    if name in filter[c] and name not in filter[w]:
                        uniqueA[name] = filter[c][name]
                    elif name in filter[w] and name not in filter[c]:
                        uniqueB[name] = filter[w][name]

                uniqueida.append(uniqueA)
                uniqueidb.append(uniqueB)

            if Q == -1:
                return ttestlist, uniqueida, uniqueidb
            if Q != -1:
                return qtestlist, uniqueida, uniqueidb

        if len(filter) == 2:

            ttestdict = {}
            qtestdict = {}
            qtestlist = []

            setstateA = set(filter[0])
            setstateB = set(filter[1])

            if Q == -1:
                for name in setstateA.intersection(setstateB):
                    FC = mean(filter[0][name][2]) / mean(filter[1][name][2])
                    ttestdict[name] = stats.ttest_ind(filter[0][name][2], filter[1][name][2])[1], filter[0][name], filter[1][name], FC
            elif Q != -1:
                for name in setstateA.intersection(setstateB):
                    ttestdict[name] = stats.ttest_ind(filter[0][name][2], filter[1][name][2])[1],
                tvals = [i[0] for i in ttestdict.values()]
                qtestlist.append(multipletests(tvals, method='fdr_bh')[1])
                g = 0
                for name in setstateA.intersection(setstateB):
                    FC = mean(filter[0][name][2]) / mean(filter[1][name][2])
                    qtestdict[name] = qtestlist[0][g], filter[0][name], filter[1][name], FC
                    g += 1

            uniqueA = {}
            uniqueB = {}

            for name in (setstateA | setstateB):
                if name in filter[0] and name not in filter[1]:
                    uniqueA[name] = filter[0][name]
                elif name in filter[1] and name not in filter[0]:
                    uniqueB[name] = filter[1][name]

            if Q == -1:
                return ttestdict, uniqueA, uniqueB
            if Q != -1:
                return qtestdict, uniqueA, uniqueB

    UpReg = []
    DownReg = []
    Unchanged = []

    def ExcelMultistateSampleOutputs(name, ttest, uniquea, uniqueb, filtered, Q, Qvalue, control, wild):

        print("Exporting Excel Files...")

        Testtype = "TTest"
        cutvalue = 0.05

        if Q == 1:
            Testtype = "Adj-QTest"
            cutvalue = Qvalue
            subdirectory = "/Adjusted Q Tests/"
            os.makedirs(DayandTime + subdirectory)

        ##DataQuality
        if Q == -1:
            subdirectory = "/DataAnalysis/"
            os.makedirs(DayandTime + subdirectory)
            for y in range(len(filtered)):
                listheadersSpC = list(str(name[y]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersNSAF = list(str(name[y]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersloge = list(str(name[y]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
                newpath = os.path.join(DayandTime + subdirectory, "DataQuality-" + str(name[y]) + ".csv")
                with open(newpath, 'w') as myfile:
                    wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                    wr.writerow(combinedheader)
                    for protein in filtered[y]:
                        SpC = list(s for s in filtered[y][protein][0])
                        NSAF = list(n for n in filtered[y][protein][2])
                        loge = list(l for l in filtered[y][protein][3])
                        mw = filtered[y][protein][1][0]
                        combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                        wr.writerow(combinedrow)
            subdirectory = "/T Tests/"
            os.makedirs(DayandTime + subdirectory)

        for x in range(len(wild)):

            UpRegbuffer = []
            DownRegbuffer = []
            Unchangedbuffer = []

            c = name.index(control)
            w = name.index(wild[x])

            ##All data
            listheadersSpC0 = list(str(name[c]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersSpC1 = list(str(name[w]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF0 = list(str(name[c]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF1 = list(str(name[w]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge0 = list(str(name[c]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge1 = list(str(name[w]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-AllData-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    SpC0 = list(s for s in ttest[x][protein][1][0])
                    SpC1 = list(s for s in ttest[x][protein][2][0])
                    NSAF0 = list(n for n in ttest[x][protein][1][2])
                    NSAF1 = list(n for n in ttest[x][protein][2][2])
                    loge0 = list(l for l in ttest[x][protein][1][4])
                    loge1 = list(l for l in ttest[x][protein][2][4])
                    mw = ttest[x][protein][1][1]
                    FC = ttest[x][protein][3]
                    tt = ttest[x][protein][0]
                    combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                    wr.writerow(combinedrow)

            ##Unchanged
            listheadersSpC0 = list(str(name[c]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersSpC1 = list(str(name[w]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF0 = list(str(name[c]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF1 = list(str(name[w]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge0 = list(str(name[c]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge1 = list(str(name[w]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Unchanged-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    if ttest[x][protein][0] > cutvalue:
                        SpC0 = list(s for s in ttest[x][protein][1][0])
                        SpC1 = list(s for s in ttest[x][protein][2][0])
                        NSAF0 = list(n for n in ttest[x][protein][1][2])
                        NSAF1 = list(n for n in ttest[x][protein][2][2])
                        loge0 = list(l for l in ttest[x][protein][1][4])
                        loge1 = list(l for l in ttest[x][protein][2][4])
                        mw = ttest[x][protein][1][1]
                        FC = ttest[x][protein][3]
                        tt = ttest[x][protein][0]
                        combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                        wr.writerow(combinedrow)
                        Unchangedbuffer.append(protein)
            Unchanged.append(Unchangedbuffer)

            ##Downregulated
            listheadersSpC0 = list(str(name[c]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersSpC1 = list(str(name[w]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF0 = list(str(name[c]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF1 = list(str(name[w]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge0 = list(str(name[c]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge1 = list(str(name[w]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Downregulated-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    if ttest[x][protein][0] < cutvalue:
                        if ttest[x][protein][3] < 1:
                            SpC0 = list(s for s in ttest[x][protein][1][0])
                            SpC1 = list(s for s in ttest[x][protein][2][0])
                            NSAF0 = list(n for n in ttest[x][protein][1][2])
                            NSAF1 = list(n for n in ttest[x][protein][2][2])
                            loge0 = list(l for l in ttest[x][protein][1][4])
                            loge1 = list(l for l in ttest[x][protein][2][4])
                            mw = ttest[x][protein][1][1]
                            FC = ttest[x][protein][3]
                            tt = ttest[x][protein][0]
                            combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                            wr.writerow(combinedrow)
                            DownRegbuffer.append(protein)
            DownReg.append(DownRegbuffer)

            ##Upregulated
            listheadersSpC0 = list(str(name[c]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersSpC1 = list(str(name[w]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF0 = list(str(name[c]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF1 = list(str(name[w]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge0 = list(str(name[c]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge1 = list(str(name[w]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Upregulated-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    if ttest[x][protein][0] < cutvalue:
                        if ttest[x][protein][3] > 1:
                            SpC0 = list(s for s in ttest[x][protein][1][0])
                            SpC1 = list(s for s in ttest[x][protein][2][0])
                            NSAF0 = list(n for n in ttest[x][protein][1][2])
                            NSAF1 = list(n for n in ttest[x][protein][2][2])
                            loge0 = list(l for l in ttest[x][protein][1][4])
                            loge1 = list(l for l in ttest[x][protein][2][4])
                            mw = ttest[x][protein][1][1]
                            FC = ttest[x][protein][3]
                            tt = ttest[x][protein][0]
                            combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                            wr.writerow(combinedrow)
                            UpRegbuffer.append(protein)
            UpReg.append(UpRegbuffer)

            ##UniqueControl
            if Q == -1:
                listheadersSpC = list(str(name[c]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersNSAF = list(str(name[c]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersloge = list(str(name[c]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
                newpath = os.path.join(DayandTime + subdirectory, "Unique-" + str(name[c]) + "-in-" + str(name[c]) + "&" + str(name[w]) + ".csv")
                with open(newpath, 'w') as myfile:
                    wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                    wr.writerow(combinedheader)
                    for protein in uniquea[x]:
                        SpC = list(s for s in uniquea[x][protein][0])
                        NSAF = list(n for n in uniquea[x][protein][2])
                        loge = list(l for l in uniquea[x][protein][3])
                        mw = uniquea[x][protein][1]
                        combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                        wr.writerow(combinedrow)

            ##UniqueWild
            if Q == -1:
                listheadersSpC = list(str(name[w]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersNSAF = list(str(name[w]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersloge = list(str(name[w]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
                newpath = os.path.join(DayandTime + subdirectory, "Unique-" + str(name[w]) +  "-in-" + str(name[w]) + "&" + str(name[c]) + ".csv")
                with open(newpath, 'w') as myfile:
                    wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                    wr.writerow(combinedheader)
                    for protein in uniqueb[x]:
                        SpC = list(s for s in uniqueb[x][protein][0])
                        NSAF = list(n for n in uniqueb[x][protein][2])
                        loge = list(l for l in uniqueb[x][protein][3])
                        mw = uniqueb[x][protein][1]
                        combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                        wr.writerow(combinedrow)

        print("Outputs Finished")
        return UpReg, DownReg, Unchanged

    def ExcelTwoSampleOutputs(name, ttest, uniquea, uniqueb, filtered, Q, Qvalue):

        print("Exporting Excel Files...")

        UpRegbuffer = []
        DownRegbuffer = []
        Unchangedbuffer = []

        Testtype = "TTest"
        cutvalue = 0.05

        if Q == 1:
            Testtype = "Adj-QTest"
            cutvalue = Qvalue
            subdirectory = "/Adjusted Q Tests/"
            os.makedirs(DayandTime + subdirectory)

        ##DataQuality
        if Q == -1:
            subdirectory = "/DataAnalysis/"
            os.makedirs(DayandTime + subdirectory)
            for y in range(len(filtered)):
                listheadersSpC = list(str(name[y]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersNSAF = list(str(name[y]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersloge = list(str(name[y]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                combinedheader =  ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
                newpath = os.path.join(DayandTime + subdirectory, "DataQuality-" + str(name[y]) + ".csv")
                with open(newpath, 'w') as myfile:
                    wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                    wr.writerow(combinedheader)
                    for protein in filtered[y]:
                        SpC = list(s for s in filtered[y][protein][0])
                        NSAF = list(n for n in filtered[y][protein][2])
                        loge = list(l for l in filtered[y][protein][3])
                        mw = filtered[y][protein][1][0]
                        combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                        wr.writerow(combinedrow)
            subdirectory = "/T Tests/"
            os.makedirs(DayandTime + subdirectory)

        ##All data
        listheadersSpC0 = list(str(name[0]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersSpC1 = list(str(name[1]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF0 = list(str(name[0]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF1 = list(str(name[1]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge0 = list(str(name[0]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge1 = list(str(name[1]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
        newpath = os.path.join(DayandTime + subdirectory, Testtype + "-AllData-" + str(name[0] + "&" + str(name[1])) + ".csv")
        with open(newpath, 'w') as myfile:
            wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
            wr.writerow(combinedheader)
            for protein in ttest:
                SpC0 = list(s for s in ttest[protein][1][0])
                SpC1 = list(s for s in ttest[protein][2][0])
                NSAF0 = list(n for n in ttest[protein][1][2])
                NSAF1 = list(n for n in ttest[protein][2][2])
                loge0 = list(l for l in ttest[protein][1][4])
                loge1 = list(l for l in ttest[protein][2][4])
                mw = ttest[protein][1][1]
                FC = ttest[protein][3]
                tt = ttest[protein][0]
                combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                wr.writerow(combinedrow)

        ##Unchanged
        listheadersSpC0 = list(str(name[0]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersSpC1 = list(str(name[1]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF0 = list(str(name[0]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF1 = list(str(name[1]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge0 = list(str(name[0]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge1 = list(str(name[1]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
        newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Unchanged-" + str(name[0] + "&" + str(name[1])) + ".csv")
        with open(newpath, 'w') as myfile:
            wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
            wr.writerow(combinedheader)
            for protein in ttest:
                if ttest[protein][0] > cutvalue:
                    SpC0 = list(s for s in ttest[protein][1][0])
                    SpC1 = list(s for s in ttest[protein][2][0])
                    NSAF0 = list(n for n in ttest[protein][1][2])
                    NSAF1 = list(n for n in ttest[protein][2][2])
                    loge0 = list(l for l in ttest[protein][1][4])
                    loge1 = list(l for l in ttest[protein][2][4])
                    mw = ttest[protein][1][1]
                    FC = ttest[protein][3]
                    tt = ttest[protein][0]
                    combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                    wr.writerow(combinedrow)
                    Unchangedbuffer.append(protein)
        Unchanged.append(Unchangedbuffer)

        ##Downregulated
        listheadersSpC0 = list(str(name[0]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersSpC1 = list(str(name[1]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF0 = list(str(name[0]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF1 = list(str(name[1]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge0 = list(str(name[0]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge1 = list(str(name[1]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
        newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Downregulated-" + str(name[0] + "&" + str(name[1])) + ".csv")
        with open(newpath, 'w') as myfile:
            wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
            wr.writerow(combinedheader)
            for protein in ttest:
                if ttest[protein][0] < cutvalue:
                    if ttest[protein][3] < 1:
                        SpC0 = list(s for s in ttest[protein][1][0])
                        SpC1 = list(s for s in ttest[protein][2][0])
                        NSAF0 = list(n for n in ttest[protein][1][2])
                        NSAF1 = list(n for n in ttest[protein][2][2])
                        loge0 = list(l for l in ttest[protein][1][4])
                        loge1 = list(l for l in ttest[protein][2][4])
                        mw = ttest[protein][1][1]
                        FC = ttest[protein][3]
                        tt = ttest[protein][0]
                        combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                        wr.writerow(combinedrow)
                        DownRegbuffer.append(protein)
        DownReg.append(DownRegbuffer)

        ##Upregulated
        listheadersSpC0 = list(str(name[0]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersSpC1 = list(str(name[1]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF0 = list(str(name[0]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF1 = list(str(name[1]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge0 = list(str(name[0]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge1 = list(str(name[1]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersNSAF0, *listheadersNSAF1, *listheadersloge0, *listheadersloge1, "Fold Change", Testtype]
        newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Upregulated-" + str(name[0] + "&" + str(name[1])) + ".csv")
        with open(newpath, 'w') as myfile:
            wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
            wr.writerow(combinedheader)
            for protein in ttest:
                if ttest[protein][0] < cutvalue:
                    if ttest[protein][3] > 1:
                        SpC0 = list(s for s in ttest[protein][1][0])
                        SpC1 = list(s for s in ttest[protein][2][0])
                        NSAF0 = list(n for n in ttest[protein][1][2])
                        NSAF1 = list(n for n in ttest[protein][2][2])
                        loge0 = list(l for l in ttest[protein][1][4])
                        loge1 = list(l for l in ttest[protein][2][4])
                        mw = ttest[protein][1][1]
                        FC = ttest[protein][3]
                        tt = ttest[protein][0]
                        combinedrow = [protein, mw, *SpC0, *SpC1, *NSAF0, *NSAF1, *loge0, *loge1, FC, tt]
                        wr.writerow(combinedrow)
                        UpRegbuffer.append(protein)
        UpReg.append(UpRegbuffer)

        ##UniqueControl
        if Q == -1:
            listheadersSpC = list(str(name[y]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF = list(str(name[y]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge = list(str(name[y]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
            newpath = os.path.join(DayandTime + subdirectory, "Unique-" + str(name[0]) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in uniquea:
                    SpC = list(s for s in uniquea[protein][0])
                    NSAF = list(n for n in uniquea[protein][2])
                    loge = list(l for l in uniquea[protein][3])
                    mw = uniquea[protein][1]
                    combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                    wr.writerow(combinedrow)

        ##UniqueWild
        if Q == -1:
            listheadersSpC = list(str(name[y]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersNSAF = list(str(name[y]) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            listheadersloge = list(str(name[y]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
            combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
            newpath = os.path.join(DayandTime + subdirectory, "Unique-" + str(name[1]) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in uniqueb:
                    SpC = list(s for s in uniqueb[protein][0])
                    NSAF = list(n for n in uniqueb[protein][2])
                    loge = list(l for l in uniqueb[protein][3])
                    mw = uniqueb[protein][1]
                    combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                    wr.writerow(combinedrow)

        return UpReg, DownReg, Unchanged

    def ExcelOneSampleOutputs(name, filtered):

        ##DataQuality

        subdirectory = "/DataAnalysis/"
        os.makedirs(DayandTime + subdirectory)

        listheadersSpC = list(str(name) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersNSAF = list(str(name) + " NSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        listheadersloge = list(str(name) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
        combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersNSAF, *listheadersloge]
        newpath = os.path.join(DayandTime + subdirectory, str(name) + ".csv")
        with open(newpath, 'w') as myfile:
            wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
            wr.writerow(combinedheader)
            for protein in filtered[0]:
                SpC = list(s for s in filtered[0][protein][0])
                NSAF = list(n for n in filtered[0][protein][2])
                loge = list(l for l in filtered[0][protein][3])
                mw = filtered[0][protein][1][0]
                combinedrow = [protein, mw, *SpC, *NSAF, *loge]
                wr.writerow(combinedrow)

    def GraphValues(filter, sets, name, unfiltered, TTest, uniquea, uniqueb):

        print("Producing graphs and figures...")

        subdirectory = "/DataAnalysis/"

        if sets == 2:
            newpath = os.path.join(DayandTime + subdirectory, str(name[0]) + " vs " + str(name[1]) + "P Value Histogram - Shared IDs.jpg")
            ttestlist = [TTest[x][0] for x in TTest.keys()]
            plt.hist(ttestlist, bins=20, histtype='stepfilled', color='b')
            plt.title("TTest for " + str(CSVInputName[0]) + " and " + str(CSVInputName[1]) + "P Value Histogram - Shared IDs")
            plt.savefig(newpath)

            plt.clf()

            newpath = os.path.join(DayandTime + subdirectory, str(name[0]) + " vs " + str(name[1]) + "P Value Histogram - Total.jpg")
            for a in range(len(uniquea) + len(uniqueb)):
                ttestlist.append(0.000444)
            plt.hist(ttestlist, bins=20, histtype='stepfilled', color='b')
            plt.title("TTest for " + str(CSVInputName[0]) + " and " + str(CSVInputName[1]) + "P Value Histogram - Total IDs")
            plt.savefig(newpath)

            plt.clf()

            newpath = os.path.join(DayandTime + subdirectory, str(name[0]) + " vs " + str(name[1]) + " Protein IDs.jpg")
            graph = venn2([set(filter[0]), set(filter[1])], (str(name[0]), str(name[1])))
            graph.get_patch_by_id('10').set_color('blue')
            graph.get_patch_by_id('01').set_color('pink')
            plt.savefig(newpath)

        for y in range(sets):
            ig, ax = plt.subplots(1, DesiredReplicateNumber, tight_layout=True, figsize=(16, 4))
            for q in range(DesiredReplicateNumber):
                values = [filter[y][x][2][q] for x in filter[y].keys()]
                N, bins, patches = ax[q].hist(values, bins=50, range=[min(values),max(values)], label=str(name[y]))
                ax[q].title.set_text("NSAF Values for " + str(name[y]) + " R" + str(q+1))
                fracs = N.astype(float) / N.max()
                norm = colors.Normalize(fracs.min(), fracs.max())
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.plasma(norm(thisfrac))
                    thispatch.set_facecolor(color)
            newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " NSAF Values for Scrappy Data.jpg")
            plt.savefig(newpath)

        for y in range(sets):
            ig, ax = plt.subplots(1, DesiredReplicateNumber, tight_layout=True, figsize=(16, 4))
            for q in range(DesiredReplicateNumber):
                values = [unfiltered[y][x][2][q] for x in unfiltered[y].keys()]
                N, bins, patches = ax[q].hist(values, bins=50, range=[min(values),max(values)], label=str(name[y]))
                ax[q].title.set_text("NSAF Values for " + str(name[y]) + " R" + str(q+1))
                fracs = N.astype(float) / N.max()
                norm = colors.Normalize(fracs.min(), fracs.max())
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.cool(norm(thisfrac))
                    thispatch.set_facecolor(color)
            newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " NSAF Values for Raw Data.jpg")
            plt.savefig(newpath)

        for y in range(sets):
            ig, ax = plt.subplots(1, DesiredReplicateNumber, tight_layout=True, figsize=(16, 4))
            for q in range(DesiredReplicateNumber):
                values = [filter[y][x][4][q] for x in filter[y].keys()]
                N, bins, patches = ax[q].hist(values, bins=100, range=[min(values),max(values)], label=str(name[y]))
                ax[q].title.set_text("log(e) Values for " + str(name[y]) + " R" + str(q+1))
                fracs = N.astype(float) / N.max()
                norm = colors.Normalize(fracs.min(), fracs.max())
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.spring(norm(thisfrac))
                    thispatch.set_facecolor(color)
            newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " Nloge Values for Scrappy Data.jpg")
            plt.savefig(newpath)
            plt.close()

    if len(CSVInputName) > 2:
        TTest, UniqueA, UniqueB = DoTTest(CSVInputName, UniqueIDsSpcNSAFFilteredDict, Control, Wild, -1)
        ExcelMultistateSampleOutputs(CSVInputName, TTest, UniqueA, UniqueB, UniqueIDsSpcNSAFFilteredDict, -1, 0, Control, Wild)

    if len(CSVInputName) == 2:
        TTest, UniqueA, UniqueB = DoTTest(UniqueIDsSpcNSAFDict, UniqueIDsSpcNSAFFilteredDict, Control, Wild, -1)
        ExcelTwoSampleOutputs(CSVInputName, TTest, UniqueA, UniqueB, UniqueIDsSpcNSAFFilteredDict, -1, 0)

    if len(CSVInputName) == 1:
        ExcelOneSampleOutputs(CSVInputName[0], UniqueIDsSpcNSAFFilteredDict)

    GraphValues(UniqueIDsSpcNSAFFilteredDict, len(CSVInputName), CSVInputName, UniqueIDsSpcNSAFDict, TTest, UniqueA, UniqueB)

    def SameSameAnalysis(a, replicates, filtered, name):

        if replicates == 6:

            subdirectory = "/SameSameAnalysis"
            os.makedirs(DayandTime + subdirectory)

            ttestlistfiltered = []
            falsediscoveryBYfiltered = []
            falsediscoveryBHfiltered = []
            bonferronifiltered = []
            Qvalue1percent = []

            for y in range(len(filtered)):

                print("Now conducting SameSame Analysis for: " + name[y])

                ttestdict = {}
                FDRBHdict = {}
                FDRBYdict = {}
                bondict = {}

                for protein in filtered[y]:
                    ttest = []
                    combination = [x for x in itertools.combinations(filtered[y][protein][2], 3)]
                    for f in range(10):
                        t = stats.ttest_ind(combination[f], combination[-(f + 1)])[1]
                        ttest.append(t)
                    averaget = mean(ttest)
                    ttestdict[protein] = (ttest, averaget)

                ttestlistfiltered.append(ttestdict)

                FDR_BY = []
                FDR_BH = []
                Bon = []

                for f in range(10):
                    fdr = {}
                    for protein in ttestdict:
                        fdr[protein] = (ttestdict[protein][0][f])
                    tvalues = [i for i in fdr.values()]
                    # FDR_BY.append(multipletests(tvalues, method='fdr_by')[1])
                    FDR_BH.append(multipletests(tvalues, method='fdr_bh')[1])
                    # Bon.append(multipletests(tvalues, method='bonferroni')[1])

                bh = []
                # by = []
                # bb = []
                for i in range(len(ttestdict)):
                    bh2 = []
                    # by2 = []
                    # bb2 = []
                    for f in range(10):
                        bh2.append(FDR_BH[f][i])
                        # by2.append(FDR_BY[f][i])
                        # bb2.append(Bon[f][i])
                    bh.append(bh2)
                    # by.append(by2)
                    # bb.append(bb2)

                i = 0
                for protein in ttestdict:
                    FDRBHdict[protein] = bh[i]
                    # FDRBYdict[protein] = by[i]
                    # bondict[protein] = bb[i]
                    i += 1

                # falsediscoveryBYfiltered.append(FDRBYdict)
                falsediscoveryBHfiltered.append(FDRBHdict)
                # bonferronifiltered.append(bondict)

                values = []
                totlen = len(ttestlistfiltered[y])

                for f in range(10):

                    valuest = []
                    # valuesfBY = []
                    valuesfBH = []
                    # valuesbon = []

                    for pqvalue in np.arange(0.05, 1.0, 0.05):

                        numbersig = len(list(1 for x in ttestlistfiltered[y].keys() if ttestlistfiltered[y][x][0][f] < pqvalue))
                        valuest.append(numbersig/totlen)

                        # numbersig2 = len(list(1 for x in falsediscoveryBYfiltered[y].keys() if falsediscoveryBYfiltered[y][x][f] < pqvalue))
                        # valuesfBY.append(numbersig2/totlen)

                        numbersig3 = len(list(1 for x in falsediscoveryBHfiltered[y].keys() if falsediscoveryBHfiltered[y][x][f] < pqvalue))
                        valuesfBH.append(numbersig3/totlen)

                        # numbersig4 = len(list(1 for x in bonferronifiltered[y].keys() if bonferronifiltered[y][x][f] < pqvalue))
                        # valuesbon.append(numbersig4/totlen)

                    # values.append([valuest] + [valuesfBY] + [valuesfBH] + [valuesbon])
                    values.append([valuest] + [valuesfBH])


                ##Pvalue Histograms
                ig, ax = plt.subplots(2,5, tight_layout=True, figsize=(18, 6))
                g = 0
                for u in range(0,2):
                    for e in range(0,5):
                        g += 1
                        histovalues = [ttestlistfiltered[y][x][0][g - 1] for x in ttestlistfiltered[y].keys()]
                        N, bins, patches = ax[u,e].hist(histovalues, bins=20, range=[min(histovalues),max(histovalues)])
                        ax[u,e].title.set_text("P Value Histogram for combination " + str(g))
                        fracs = N.astype(float) / N.max()
                        norm = colors.Normalize(fracs.min(), fracs.max())
                        for thisfrac, thispatch in zip(fracs, patches):
                            color = plt.cm.winter(norm(thisfrac))
                            thispatch.set_facecolor(color)
                newpath = os.path.join(DayandTime + subdirectory, "Ten-way Scrappy Comparison for " + str(name[y]) + ".jpg")
                plt.savefig(newpath)

                #Ttests FDR
                testnames = ["Ttest", "BH"]
                for testtype in range(2):
                    ig, ax = plt.subplots(2,5, tight_layout=True, figsize=(18, 6))
                    categories = np.arange(0.05, 1, 0.05)
                    g = 0
                    bar_width = 0.03
                    for u in range(0,2):
                        for e in range(0,5):
                            g += 1
                            fdrvalues = [x*100 for x in values[g-1][testtype]]
                            ax[u,e].bar(categories, fdrvalues, bar_width)
                            ax[u,e].title.set_text(str(testnames[testtype]) + " FDR" + str(g))
                    newpath = os.path.join(DayandTime + subdirectory, "Ten-way " + str(testnames[testtype]) + " FDR Calc for " + str(name[y]) + ".jpg")
                    plt.savefig(newpath)

                #Q value average from BH test
                ig, ax = plt.subplots(tight_layout=True)
                categories = np.arange(0.05, 1, 0.05)
                bar_width = 0.03
                fdrvalues = []
                chartvalues = []
                for u in range(10):
                    fdrvalues.append([x * 100 for x in values[u][1]])
                for x in range(len(fdrvalues[0])):
                    chartvalues.append(mean(y[x] for y in fdrvalues))
                ax.bar(categories, chartvalues, bar_width)
                ax.title.set_text("Average Q value FDR")
                newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " Average Q value FDR.jpg")
                plt.savefig(newpath)

                Qvalue1percent.append(bisect.bisect_left(chartvalues, 1.0)*0.05)

            return mean(Qvalue1percent)

    def PvsQGraph(name, UpReg, DownReg, Unchanged, control, wild):

        subdirectory = "/PvQ-Histograms"
        os.makedirs(DayandTime + subdirectory)

        print("Producing P vs Q value venn diagrams")

        if len(name) >= 3:

            i = 0

            for x in range(len(wild)):

                w = int((len(UpReg) / 2) + i)
                i += 1

                plt.clf()
                ##Up
                newpath = os.path.join(DayandTime + subdirectory, "P vs Q Upregulated" + str(control) + " vs " + str(wild[x]) + ".jpg")
                graph = venn2([set(UpReg[x]), set(UpReg[w])], ("P", "Q"))
                graph.get_patch_by_id('10').set_color('blue')
                graph.get_patch_by_id('01').set_color('pink')
                plt.savefig(newpath)

                plt.clf()
                # Down
                newpath = os.path.join(DayandTime + subdirectory, "P vs Q Downregulated" + str(control) + " vs " + str(wild[x]) + ".jpg")
                graph = venn2([set(DownReg[x]), set(DownReg[w])], ("P", "Q"))
                graph.get_patch_by_id('10').set_color('blue')
                graph.get_patch_by_id('01').set_color('pink')
                plt.savefig(newpath)

                plt.clf()
                # Unchanged
                newpath = os.path.join(DayandTime + subdirectory, "P vs Q Unchanged" + str(control) + " vs " + str(wild[x]) + ".jpg")
                graph = venn2([set(Unchanged[x]), set(Unchanged[w])], ("P", "Q"))
                graph.get_patch_by_id('10').set_color('yellow')
                graph.get_patch_by_id('01').set_color('grey')
                plt.savefig(newpath)

        if len(name) < 3:

            plt.clf()
            ##Up
            newpath = os.path.join(DayandTime + subdirectory, "P vs Q Upregulated.jpg")
            graph = venn2([set(UpReg[0]), set(UpReg[1])], ("P", "Q"))
            graph.get_patch_by_id('10').set_color('blue')
            graph.get_patch_by_id('01').set_color('pink')
            plt.savefig(newpath)

            plt.clf()
            #Down
            newpath = os.path.join(DayandTime + subdirectory, "P vs Q Downregulated.jpg")
            graph = venn2([set(DownReg[0]), set(DownReg[1])], ("P", "Q"))
            graph.get_patch_by_id('10').set_color('blue')
            graph.get_patch_by_id('01').set_color('pink')
            plt.savefig(newpath)

            plt.clf()
            #Unchanged
            newpath = os.path.join(DayandTime + subdirectory, "P vs Q Unchanged.jpg")
            graph = venn2([set(Unchanged[0]), set(Unchanged[1])], ("P", "Q"))
            graph.get_patch_by_id('10').set_color('yellow')
            graph.get_patch_by_id('01').set_color('grey')
            plt.savefig(newpath)

    if DesiredReplicateNumber == 6:
        QValue = SameSameAnalysis(UniqueIDsSpcNSAFDict, DesiredReplicateNumber, UniqueIDsSpcNSAFFilteredDict, CSVInputName)
        if len(CSVInputName) == 2:
            QTest, UniqueA, UniqueB = DoTTest(CSVInputName, UniqueIDsSpcNSAFFilteredDict, Control, Wild, QValue)
            UpReg, DownReg, Unchanged = ExcelTwoSampleOutputs(CSVInputName, QTest, UniqueA, UniqueB, UniqueIDsSpcNSAFFilteredDict, 1, QValue)
            PvsQGraph(CSVInputName, UpReg, DownReg, Unchanged, 0, 0)
        if len(CSVInputName) > 2:
            QTest, UniqueA, UniqueB = DoTTest(CSVInputName, UniqueIDsSpcNSAFFilteredDict, Control, Wild, QValue)
            UpReg, DownReg, Unchanged = ExcelMultistateSampleOutputs(CSVInputName, QTest, UniqueA, UniqueB, UniqueIDsSpcNSAFFilteredDict, 1, QValue, Control, Wild)
            PvsQGraph(CSVInputName, UpReg, DownReg, Unchanged, Control, Wild)

    endtime = time.perf_counter()
    print("Total Processing time: " + str(endtime - starttime))

    shutil.make_archive(os.path.dirname(os.path.realpath(__file__)) + "/Output", 'zip', DayandTime)
    shutil.rmtree(DayandTime)

if __name__ == "__main__":
    print("Python Version: " + str(sys.version))
    main()
