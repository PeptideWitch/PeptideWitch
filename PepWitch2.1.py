import numpy as np
import pandas as pd
import glob
import os.path
import time
import datetime
import itertools
import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib_venn import venn2
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests
from statistics import mean
import bisect
from sys import argv
import shutil
import sys
from sklearn.decomposition import PCA as sklearnPCA
from sklearn import preprocessing
from sklearn.preprocessing import normalize

def main(mode):

    #Setting up a few housekeeping items.
    def Housekeeper():
        np.set_printoptions(precision=2)
        tTest = False
        uniqueA = False
        uniqueB = False
        upReg = []
        downReg = []
        unchanged = []
        return tTest, uniqueA, uniqueB, upReg, downReg, unchanged
    TTest, UniqueA, UniqueB, UpReg, DownReg, Unchanged = Housekeeper()

    #A function that sets up all the info needed for PepWitch to collect user data names, etc.
    def GetInputName():
        filenames = 0
        control = 0
        wild = 0
        replicates = []
        dis = 0
        minspc = 5
        specfrac = 0.1
        pval = 0.05
        engine = "No"

        print("############")
        print("Welcome to PyWitch!")
        print("You can either import a single state for analysis or perform a multi-state comparison.")
        print("For a two-state comparison, please enter the control state first, then the wild type.")
        print("############")

        while filenames == 0:
            print("Provide the name of the file stems (the text before the -R1.csv tail) that you want to import, separated by a space.")
            print("Example: Control Treated1 Treated2 Treated3.")
            filenames = input("Make sure that your files end in -R#.csv, but don't type that here: ").split()

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
                                            control = int(controlstate)-1
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

        defaultvals = input("[Y/N] Do you want to use the default value set for Disregard(0)?").upper()
        if defaultvals == "N":
            while dis != 0 or 1:
                dis = int(input("Disregard (0 or 1): "))
        while engine != "GPM" or "PD":
            engine = input("Engine ~ GPM or PD: ")
            if engine == "GPM" or "PD":
                break

        return filenames, control, wild, replicatecarry, dis, minspc, specfrac, pval, engine

    #Checks to see whether the user is coming through peptidewitch.com or through a local client
    if mode == "online":
        print("We're online. You shouldn't see this, but the log will record it.")
        WorkingDirectory = argv[6]
        matplotlib.use("Agg")
        Disregard = int(argv[2])
        MinSpc = float(argv[1])
        SpectralFraction = float(argv[5])
        SameSamePValue = float(argv[3])
        Engine = argv[7]
        CSVInputName = argv[8:]
        Control = CSVInputName[0]
        if len(CSVInputName) > 1:
            Wild = CSVInputName[1:]
        else:
            Wild = CSVInputName[1]
        DesiredReplicateNumber = int(argv[4])
    else:
        WorkingDirectory = os.path.dirname(os.path.realpath(__file__))
        CSVInputName, Control, Wild, DesiredReplicateNumber, Disregard, MinSpc, SpectralFraction, SameSamePValue, Engine = GetInputName()

    #Establishes directory paths and begins timer
    DayandTime = os.path.dirname(os.path.realpath(__file__)) + "/output/" + str(datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S')) + '-' + '-'.join(CSVInputName) + "-" + str(MinSpc)
    os.makedirs(DayandTime)
    starttime = time.perf_counter()

    ## This class imports the .csv files into our workflow and makes a list of the replicates.
    def ImportProteinFiles(inputnames, engine):
        print("Importing .csv files...")

        filenames = inputnames
        TotalDataIDsL = []
        UniqueNames = []
        TotalDataIDsD = []

        for name in filenames:

            if engine == "GPM":
                indexes = ['Identifier', 'rI', 'Mr (kDa)', 'log(e)']
                try:
                    headers = pd.read_csv(str(WorkingDirectory) + "/" + str(name) + "-R1.csv", nrows=1).columns.tolist()
                    total = [pd.read_csv(x, delimiter = ',', header = 'infer',  usecols = [int(headers.index(indexes[0])),int(headers.index(indexes[1])),int(headers.index(indexes[2])),int(headers.index(indexes[3]))], quotechar = '"' , dtype = {'Identifier':'U','rl':'float32','Mr (kDa)':'float32','log(e)':'float32'}).fillna(0.01).values.tolist()  for x in glob.glob(str(WorkingDirectory) + "/" + str(name) + "*.csv")] # Imports Protein Name, Spectral Count, Mol. Weight, score
                except ValueError:
                    print("######################################")
                    print("Error #4: Your CSV file is empty or incompatible.")
                    exit()

            if engine == "PD":
                indexes = ["Accession", "# Peptides", "MW [kDa]", "Score"]
                try:
                    headers = pd.read_csv(str(WorkingDirectory) + "/" + str(name) + "-R1.csv", nrows=1).columns.tolist()
                    total = [pd.read_csv(x, delimiter = ',', header = 'infer', usecols = [int(headers.index(indexes[0])),int(headers.index(indexes[1])),int(headers.index(indexes[2])),int(headers.index(indexes[3]))], quotechar = '"' , dtype = {'Accession':'U','# Peptides':'float32','MW [kDa]':'float32','Score':'float32'}).fillna(0.01).values.tolist() for x in glob.glob(str(WorkingDirectory) + "/" + str(name) + "*.csv")] # Imports Protein Name, Peptides, Mol. Weight, score
                except ValueError:
                    print("######################################")
                    print("Error #4: Your CSV file is empty or incompatible.")
                    exit()

            totaldict = [dict((row[0], row[1:]) for row in t) for t in total]
            uniquenameset = set(x for data in totaldict for x in data)  # (for data in totaldict: for x in data: x)

            TotalDataIDsL.append(total)
            TotalDataIDsD.append(totaldict)
            UniqueNames.append(uniquenameset)

        return TotalDataIDsL, UniqueNames, TotalDataIDsD
    TotalDataIDsList, UniqueNameList, TotalDataIDsDictionary, = ImportProteinFiles(CSVInputName, Engine)

    ## This will take the Total ID data and produce NSAF values based off the original state
    def CalculatelnNSAF(flatlist, spec):
        print("Calculating lnNSAF")

        newcopy = flatlist

        Collectivemolw = []
        Collectivespc = []
        TotalDict = []
        ZeroValueMolW = []
        ZeroSAF = []
        ZerolnNSAF = []

        for state in newcopy: #For control, wild, etc

            collectivemolwstate = []
            collectivespcstate = []
            collectiveSAFstate = []
            replicatesumSAF = []
            stateDict = []
            zerovaluemolwstate = []
            zeroSAFstate = []
            zerolnNSAFstate = []

            # Here, we collect the total MolW and SpC counts from the raw data.
            for replicate in range(len(state)):
                replicatesummolw = 0
                replicatesumspc = 0
                averagelist = []
                for protein in state[replicate]:
                    molw = (protein[3] * 1000) / 110  # Adds all molecular weights into aa
                    replicatesummolw += molw
                    averagelist.append(molw)
                    replicatesumspc += (protein[2])  # Adds together all the SpCs
                collectivemolwstate.append(replicatesummolw)
                collectivespcstate.append(replicatesumspc)
                zerovaluemolwstate.append(np.median(averagelist))
            Collectivemolw.append(collectivemolwstate)
            Collectivespc.append(collectivespcstate)
            ZeroValueMolW.append(zerovaluemolwstate)

            # SAF calculation
            for replicate in range(len(state)):
                replicatesumSAFbuffer = 0
                for index, protein in enumerate(state[replicate]):
                    bufferval = ( (protein[2] / collectivespcstate[replicate]) / ( (1000 * protein[3] / 110) / (collectivemolwstate[replicate])) ) # Here, the SAF is calculated with relation to the individual replicate snapshot [CollectMolW/SpC], not the filtered snapshot
                    state[replicate][index].append(bufferval)
                    replicatesumSAFbuffer += bufferval
                zeroSAFstate.append(spec / zerovaluemolwstate[replicate])
                replicatesumSAF.append(replicatesumSAFbuffer)
            collectiveSAFstate.append(replicatesumSAF)
            ZeroSAF.append(zeroSAFstate)

            # lnNSAF calculation
            for replicate in range(len(state)):
                zerolnNSAFstate.append(np.log(zeroSAFstate[replicate] / replicatesumSAF[replicate]))
                for index, protein in enumerate(state[replicate]):
                    bufferval = np.log(protein[4] / replicatesumSAF[replicate])
                    state[replicate][index].append(bufferval)
            ZerolnNSAF.append(zerolnNSAFstate)

            #Converting object to dictionary
            for replicate in range(len(state)):
                TotalDictbuffer = {}
                for protein in state[replicate]:
                    TotalDictbuffer[protein[0]] = list(protein[1:])
                stateDict.append(TotalDictbuffer)
            TotalDict.append(stateDict)

        return TotalDict, ZeroValueMolW, ZeroSAF, ZerolnNSAF
    AllDataDict, ZeroMWL, ZeroSAFL, ZerolnNSAFL = CalculatelnNSAF(TotalDataIDsList, SpectralFraction)

    # This will take the NSAF data and, using the Unique values from all replicates, produce a filtered Scrappy Rules Dictionary.
    def CollectUniqueIdsAndSpectra(unique, alldata, spec, engine, disregard, zeromw, zerosaf, zerolnnsaf):
        print("Finding Unique IDs...")

        TotalFilter = []

        for y in range(len(unique)):  # for every test state [control, wild]
            dummy = [spec, spec, spec, zerosaf[y], zerolnnsaf[y]]
            filteredvaluesdict = {}  # This is Scrappy, basically. Conducts as minimum spectral count and disregard replicate drop function
            for name in unique[y]:
                if ([x[1] for x in (data.get(name, dummy) for data in alldata[y])].count(spec) <= disregard) and (sum(x[1] for x in (data.get(name, dummy) for data in alldata[y])) >= MinSpc):
                    filteredvaluesdict[name] = [
                                [x[1] if x[1] != 0 else spec for x in (replicate.get(name, dummy) for replicate in alldata[y])], # SpC
                                [b for b in (next(x for x in (replicate.get(name, [spec])[2:3] for replicate in alldata[y]) if x))], # Molw
                                [i[0] if i[0] != 0 else spec for i in (replicate.get(name, dummy) for replicate in alldata[y])], # Score
                                [x[3][i] if type(x[3]) == list else x[3] for i,x in ((index, replicate.get(name, dummy)) for index, replicate in enumerate(alldata[y]))], # SAF
                                [x[4][i] if type(x[4]) == list else x[4] for i,x in ((index, replicate.get(name, dummy)) for index, replicate in enumerate(alldata[y]))]  # lnNSAF
                                                ]
            TotalFilter.append(filteredvaluesdict)

        return TotalFilter

    ScrappyIDsDict = CollectUniqueIdsAndSpectra(UniqueNameList, AllDataDict, SpectralFraction, Engine, Disregard, ZeroMWL, ZeroSAFL, ZerolnNSAFL)

    def DoTTest(filter, control, wild, Q, inputnames):

        if Q == -1:
            print("Now conducting TTests for all proteins...")
        elif Q != -1:
            print("Now conducting Adjusted QTests for all proteins...")

        if len(filter) > 2:

            uniqueida = []
            uniqueidb = []
            ttestlist = []
            qtestlist = []

            for state in wild: #For a multistate analysis

                qtestlista = []
                ttestdict = {}
                qtestdict = {}

                c = control
                w = inputnames.index(state)

                setstateControl = set(filter[c])
                setstateWild = set(filter[w])

                if Q == -1:
                    for name in setstateControl.intersection(setstateWild):
                        FC = mean(filter[w][name][3]) / mean(filter[c][name][3])
                        ttestdict[name] = stats.ttest_ind(filter[c][name][4], filter[w][name][4])[1], filter[c][name], filter[w][name], FC
                    ttestlist.append(ttestdict)
                elif Q != -1:
                    for name in setstateControl.intersection(setstateWild):
                        ttestdict[name] = stats.ttest_ind(filter[c][name][4], filter[w][name][4])[1],
                    tvals = [i[0] for i in ttestdict.values()]
                    qtestlista.append(multipletests(tvals, method='fdr_bh')[1])
                    g = 0
                    for name in setstateControl.intersection(setstateWild):
                        FC = mean(filter[w][name][3]) / mean(filter[c][name][3])
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
                    FC = mean(filter[1][name][3]) / mean(filter[0][name][3]) #Fold change calculated from SAF, not lnNSAF
                    ttestdict[name] = stats.ttest_ind(filter[0][name][4], filter[1][name][4])[1], filter[0][name], filter[1][name], FC  #TTest calculated from lnNSAF
            elif Q != -1:
                for name in setstateA.intersection(setstateB):
                    ttestdict[name] = stats.ttest_ind(filter[0][name][4], filter[1][name][4])[1],
                tvals = [i[0] for i in ttestdict.values()]
                qtestlist.append(multipletests(tvals, method='fdr_bh')[1])
                g = 0
                for name in setstateA.intersection(setstateB):
                    FC = mean(filter[1][name][3]) / mean(filter[0][name][3]) #Fold change calculated from SAF, not lnNSAF
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

    def PrincipleComponent(inputnames, scrappy, mode, con, wil, reps):

        safoptions = [3,4]
        saftext = ["SAF", "lnNSAF"]
        typetext = ["Replicates", "Proteins"]

        if mode == "Data":
            print("Now conducting PCA for Individual Replicates")
            for index, option in enumerate(safoptions):
                for state in range(len(inputnames)):
                    arr = []
                    IDs = list(str(i) for i in scrappy[state])
                    for protein in scrappy[state]:
                        arr.append(scrappy[state][protein][option]) #SAF data first, then lnNSAF second, filtering out impute values

                    for repeat in range(0,2):
                        if repeat == 0:
                            data = preprocessing.StandardScaler().fit_transform(arr).transpose() #.transpose = Replicate Analysis
                        if repeat == 1:
                            data = preprocessing.StandardScaler().fit_transform(arr) # non transpose = Protein Analysis
                        data2 = sklearnPCA(reps).fit_transform(data)

                        # if sklearnPCA(2).fit(data).explained_variance_ratio_[0] > .5:
                        #     abspc1 = list(abs(i[0]) for i in data2)
                        #     outliers = list(IDs[index] for index, item in enumerate(abspc1) if stats.percentileofscore(abspc1, item) > 95)

                        principalDf =  pd.DataFrame(data = data2[:,0:2], columns=['PC1', 'PC2'])
                        fig = plt.figure(figsize = (6,6))
                        ax = fig.add_subplot(1,1,1)
                        ax.set_xlabel("PC1", fontsize=14)
                        ax.set_ylabel("PC2", fontsize=14)
                        ax.text(-1,1.5, str(sklearnPCA(reps).fit(data).explained_variance_ratio_), fontsize=7)
                        ax.set_title("2 Component PCA for " + str(inputnames[state] + " " + str(saftext[index]) + " " + str(typetext[repeat])), fontsize=14)
                        ax.scatter(principalDf['PC1'], principalDf['PC2'], s=50)
                        ax.axhline(0,0,1, c = 'Green', linewidth=3)
                        ax.axvline(0,0,1, c = 'Purple',linewidth=3)
                        ax.grid()
                        subdirectory = "/DataAnalysis/"
                        newpath = os.path.join(DayandTime + subdirectory, str(inputnames[state] + " " + str(saftext[index]) + " " + str(typetext[repeat]) + " PCA.jpg"))
                        plt.savefig(newpath)
                        plt.clf()

        if mode == "TTest":

            if len(inputnames) == 2:
                scrappynew = [scrappy, scrappy]
            if len(inputnames) > 2:
                scrappynew = scrappy

            if len(inputnames) == 2:
                wild = [wil, wil]
            if len(inputnames) > 2:
                wild = wil

            print("Now conducting PCA for Ttest states")
            for state in range(len(inputnames)-1):
                arr = []
                IDs = list(str(i) for i in scrappynew[state])
                for protein in scrappynew[state]:
                    arr.append(scrappynew[state][protein][1][4] + scrappynew[state][protein][2][4])  # lnNSAF data

                data = preprocessing.StandardScaler().fit_transform(arr).transpose()
                data2 = sklearnPCA(reps).fit_transform(data)

                if sklearnPCA(2).fit(data).explained_variance_ratio_[0] > .5:
                    abspc1 = list(abs(i[0]) for i in data2)
                    outliers = list(
                        IDs[index] for index, item in enumerate(abspc1) if stats.percentileofscore(abspc1, item) > 95)

                principalDf = pd.DataFrame(data=data2[:,0:2], columns=['PC1', 'PC2'])
                targetDf = pd.DataFrame(np.array( list(inputnames[con] for item in range(reps)) + list(wild[state] for item in range(reps)) ), columns=["State"])
                finalDf = pd.concat([principalDf, targetDf], axis=1)
                fig = plt.figure(figsize=(6, 6))
                ax = fig.add_subplot(1, 1, 1)
                ax.set_xlabel("PC1", fontsize=14)
                ax.set_ylabel("PC2", fontsize=14)
                ax.text(-1, 1.5, str(sklearnPCA(reps).fit(data).explained_variance_ratio_))
                ax.set_title("2 Component PCA Fit Transform for TTest " + str(inputnames[con]) + " vs " + str(wild[state]), fontsize=12)
                targets = [inputnames[con], wild[state]]
                colours = ['r', 'b']
                for target, colour in zip(targets, colours):
                    indicesToKeep = finalDf['State'] == target
                    ax.scatter(finalDf.loc[indicesToKeep, 'PC1'], finalDf.loc[indicesToKeep, 'PC2'], c=colour, s=40)
                ax.legend(targets)
                ax.axhline(0, 0, 1, c='Green', linewidth=3)
                ax.axvline(0, 0, 1, c='Purple', linewidth=3)
                ax.grid()
                subdirectory = "/T Tests/"
                newpath = os.path.join(DayandTime + subdirectory, str(inputnames[con]) + " vs " + str(wild[state]) + " PCA.jpg")
                plt.savefig(newpath)
                plt.clf()

    def Volcano(inputnames, ttest, con, wil):

        print("Making a Volcano plot")

        if len(inputnames) == 2:
            wild = [wil, wil]
        if len(inputnames) > 2:
            wild = wil

        if len(inputnames) == 2:
            ttestnew = [ttest, ttest]
        if len(inputnames) > 2:
            ttestnew = ttest

        for state in range(len(inputnames)-1):

            dict = {"ID": list(i for i in ttestnew[state]), "P-Value": list(ttestnew[state][i][0] for i in ttestnew[state]), "Fold Change": list(ttestnew[state][i][3] for i in ttestnew[state])}
            dict2 = {"P-Value": dict["P-Value"], "Fold Change": dict["Fold Change"]}
            df1 = pd.DataFrame(dict2["P-Value"], columns=["P-Value"])
            df2 = pd.DataFrame(dict2["Fold Change"], columns=["Fold Change"])
            transform_df1 = pd.DataFrame(abs(np.log10(df1["P-Value"])), columns=["P-Value"])
            transform_df2 = pd.DataFrame(np.log2(df2["Fold Change"]), columns=["Fold Change"])
            transform_dffinal = pd.DataFrame()
            transform_dffinal["P-Value"] = transform_df1
            transform_dffinal["Fold Change"] = transform_df2

            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel("log2 FC", fontsize=14)
            ax.set_ylabel("abs(log10) P Value", fontsize=14)
            ax.set_title("Volcano Plot for " + str(inputnames[con]) + " vs " + str(wild[state]), fontsize=14)
            ax.scatter(transform_dffinal["Fold Change"], transform_dffinal["P-Value"], s=2)
            ax.grid()

            subdirectory = "/T Tests/"
            newpath = os.path.join(DayandTime + subdirectory, str(inputnames[con]) + " vs " + str(wild[state]) + " Volcano.jpg")
            plt.savefig(newpath)
            plt.clf()

    def HeatMap(inputnames, ttest):

        print("Making a HeatMap")

        if ttest is not list:
            ttestnew = [ttest, ttest]
        else:
            ttestnew = ttest

        for state in range(len(inputnames)-1):

            heatdataControl = []
            heatdataWild = []
            heatdata = []

            for protein in ttestnew[state]:
                heatdataControl.append(ttestnew[state][protein][1][4]) # State 1 SAF
                heatdataWild.append(ttestnew[state][protein][2][4]) # State 2 SAF

            for row in range(len(heatdataWild)):
                heatdata.append(heatdataControl[row] + heatdataWild[row])

            fig = plt.figure(figsize=(60, 60))
            ax = fig.add_subplot(1, 1, 1)
            ax.imshow(heatdata)
            ax.set_title("HeatMap for Control and Wild")
            ax.set_aspect('auto')
            plt.show()
            plt.clf()

    def ExcelSampleOutputs(name, tteststart, uniqueastart, uniquebstart, Q, Qvalue, control, wild, Samecut, repnum):

        print("Exporting Two Sample A vs B Files...")

        if Q == 1:
            Testtype = "Adj-QTest"
            cutvalue = Qvalue
            subdirectory = "/Adjusted Q Tests/"
            os.makedirs(DayandTime + subdirectory)

        elif Q == -1:
            Testtype = "TTest"
            cutvalue = Samecut
            subdirectory = "/T Tests/"
            os.makedirs(DayandTime + subdirectory)

        for x in range(len(name)-1): # So if there are 3 states, the ttests will be conducted in two rounds

            UpRegbuffer = []
            DownRegbuffer = []
            Unchangedbuffer = []

            #This piece of code here helps the test parse. Ttest outputs will either be in dictionaries or lists depending
            #on the number of states input. So here, we make a dummy list for the A vs B, thereby allowing the code beneath
            #to run in a general sense.
            if len(name) > 2:
                w = name.index(wild[x])
                c = control # An integer
                ttest = tteststart
                uniquea = uniqueastart
                uniqueb = uniquebstart
            else:
                w = 1
                c = 0
                ttest = [tteststart, tteststart]
                uniquea = [uniqueastart, uniqueastart]
                uniqueb = [uniquebstart, uniquebstart]

            ##All data
            listheadersSpC0 = list(str(name[c]) + " SpC Control R" + str(x + 1) for x in range(0, repnum))
            listheadersSpC1 = list(str(name[w]) + " SpC Treatment R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF0 = list(str(name[c]) + " NSAF Control R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF1 = list(str(name[w]) + " NSAF Treatment R" + str(x + 1) for x in range(0, repnum))
            listheadersloge0 = list(str(name[c]) + " -log(e) Control R" + str(x + 1) for x in range(0, repnum))
            listheadersloge1 = list(str(name[w]) + " -log(e) Treatment R" + str(x + 1) for x in range(0, repnum))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersloge0, *listheadersloge1, *listheaderslnNSAF0, *listheaderslnNSAF1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-AllData-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    SpC0 = list(s for s in ttest[x][protein][1][0])
                    SpC1 = list(s for s in ttest[x][protein][2][0])
                    NSAF0 = list(n for n in ttest[x][protein][1][4])
                    NSAF1 = list(n for n in ttest[x][protein][2][4])
                    loge0 = list(l for l in ttest[x][protein][1][2])
                    loge1 = list(l for l in ttest[x][protein][2][2])
                    mw = ttest[x][protein][1][1][0]
                    FC = ttest[x][protein][3]
                    tt = ttest[x][protein][0]
                    combinedrow = [protein, mw, *SpC0, *SpC1, *loge0, *loge1, *NSAF0, *NSAF1, FC, tt]
                    wr.writerow(combinedrow)

            ##Unchanged
            listheadersSpC0 = list(str(name[c]) + " SpC Control R" + str(x + 1) for x in range(0, repnum))
            listheadersSpC1 = list(str(name[w]) + " SpC Treatment R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF0 = list(str(name[c]) + " NSAF Control R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF1 = list(str(name[w]) + " NSAF Treatment R" + str(x + 1) for x in range(0, repnum))
            listheadersloge0 = list(str(name[c]) + " -log(e) Control R" + str(x + 1) for x in range(0, repnum))
            listheadersloge1 = list(str(name[w]) + " -log(e) Treatment R" + str(x + 1) for x in range(0, repnum))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersloge0, *listheadersloge1, *listheaderslnNSAF0, *listheaderslnNSAF1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Unchanged-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    if ttest[x][protein][0] > cutvalue:
                        SpC0 = list(s for s in ttest[x][protein][1][0])
                        SpC1 = list(s for s in ttest[x][protein][2][0])
                        NSAF0 = list(n for n in ttest[x][protein][1][4])
                        NSAF1 = list(n for n in ttest[x][protein][2][4])
                        loge0 = list(l for l in ttest[x][protein][1][2])
                        loge1 = list(l for l in ttest[x][protein][2][2])
                        mw = ttest[x][protein][1][1][0]
                        FC = ttest[x][protein][3]
                        tt = ttest[x][protein][0]
                        combinedrow = [protein, mw, *SpC0, *SpC1, *loge0, *loge1, *NSAF0, *NSAF1, FC, tt]
                        wr.writerow(combinedrow)
                        Unchangedbuffer.append(protein)
            Unchanged.append(Unchangedbuffer)

            ##Downregulated
            listheadersSpC0 = list(str(name[c]) + " SpC Control R" + str(x + 1) for x in range(0, repnum))
            listheadersSpC1 = list(str(name[w]) + " SpC Treatment R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF0 = list(str(name[c]) + " NSAF Control R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF1 = list(str(name[w]) + " NSAF Treatment R" + str(x + 1) for x in range(0, repnum))
            listheadersloge0 = list(str(name[c]) + " -log(e) Control R" + str(x + 1) for x in range(0, repnum))
            listheadersloge1 = list(str(name[w]) + " -log(e) Treatment R" + str(x + 1) for x in range(0, repnum))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersloge0, *listheadersloge1, *listheaderslnNSAF0, *listheaderslnNSAF1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Downregulated-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    if ttest[x][protein][0] < cutvalue:
                        if ttest[x][protein][3] < 1:
                            SpC0 = list(s for s in ttest[x][protein][1][0])
                            SpC1 = list(s for s in ttest[x][protein][2][0])
                            NSAF0 = list(n for n in ttest[x][protein][1][4])
                            NSAF1 = list(n for n in ttest[x][protein][2][4])
                            loge0 = list(l for l in ttest[x][protein][1][2])
                            loge1 = list(l for l in ttest[x][protein][2][2])
                            mw = ttest[x][protein][1][1][0]
                            FC = ttest[x][protein][3]
                            tt = ttest[x][protein][0]
                            combinedrow = [protein, mw, *SpC0, *SpC1, *loge0, *loge1, *NSAF0, *NSAF1, FC, tt]
                            wr.writerow(combinedrow)
                            DownRegbuffer.append(protein)
            DownReg.append(DownRegbuffer)

            ##Upregulated
            listheadersSpC0 = list(str(name[c]) + " SpC Control R" + str(x + 1) for x in range(0, repnum))
            listheadersSpC1 = list(str(name[w]) + " SpC Treatment R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF0 = list(str(name[c]) + " NSAF Control R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF1 = list(str(name[w]) + " NSAF Treatment R" + str(x + 1) for x in range(0, repnum))
            listheadersloge0 = list(str(name[c]) + " -log(e) Control R" + str(x + 1) for x in range(0, repnum))
            listheadersloge1 = list(str(name[w]) + " -log(e) Treatment R" + str(x + 1) for x in range(0, repnum))
            combinedheader = ["Protein", "MW", *listheadersSpC0, *listheadersSpC1, *listheadersloge0, *listheadersloge1, *listheaderslnNSAF0, *listheaderslnNSAF1, "Fold Change", Testtype]
            newpath = os.path.join(DayandTime + subdirectory, Testtype + "-Upregulated-" + str(name[c] + "&" + str(name[w])) + ".csv")
            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in ttest[x]:
                    if ttest[x][protein][0] < cutvalue:
                        if ttest[x][protein][3] > 1:
                            SpC0 = list(s for s in ttest[x][protein][1][0])
                            SpC1 = list(s for s in ttest[x][protein][2][0])
                            NSAF0 = list(n for n in ttest[x][protein][1][4])
                            NSAF1 = list(n for n in ttest[x][protein][2][4])
                            loge0 = list(l for l in ttest[x][protein][1][2])
                            loge1 = list(l for l in ttest[x][protein][2][2])
                            mw = ttest[x][protein][1][1][0]
                            FC = ttest[x][protein][3]
                            tt = ttest[x][protein][0]
                            combinedrow = [protein, mw, *SpC0, *SpC1, *loge0, *loge1, *NSAF0, *NSAF1, FC, tt]
                            wr.writerow(combinedrow)
                            UpRegbuffer.append(protein)
            UpReg.append(UpRegbuffer)

            ##UniqueControl
            if Q == -1:
                listheadersSpC = list(str(name[c]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheaderslnNSAF = list(str(name[c]) + " lnNSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersloge = list(str(name[c]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersloge, *listheaderslnNSAF, ]
                newpath = os.path.join(DayandTime + subdirectory, "Unique-" + str(name[c]) + ".csv")
                with open(newpath, 'w') as myfile:
                    wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                    wr.writerow(combinedheader)
                    for protein in uniquea[x]:
                        SpC = list(s for s in uniquea[x][protein][0])
                        NSAF = list(n for n in uniquea[x][protein][4])
                        loge = list(l for l in uniquea[x][protein][2])
                        mw = uniquea[x][protein][1][0]
                        combinedrow = [protein, mw, *SpC, *loge, *NSAF]
                        wr.writerow(combinedrow)

            ##UniqueWild
            if Q == -1:
                listheadersSpC = list(str(name[w]) + " SpC R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheaderslnNSAF = list(str(name[w]) + " lnNSAF R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                listheadersloge = list(str(name[w]) + " -log(e) R" + str(x + 1) for x in range(0, DesiredReplicateNumber))
                combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersloge, *listheaderslnNSAF]
                newpath = os.path.join(DayandTime + subdirectory, "Unique-" + str(name[w]) + ".csv")
                with open(newpath, 'w') as myfile:
                    wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                    wr.writerow(combinedheader)
                    for protein in uniqueb[x]:
                        SpC = list(s for s in uniqueb[x][protein][0])
                        NSAF = list(n for n in uniqueb[x][protein][4])
                        loge = list(l for l in uniqueb[x][protein][2])
                        mw = uniqueb[x][protein][1][0]
                        combinedrow = [protein, mw, *SpC, *loge, *NSAF]
                        wr.writerow(combinedrow)

        return UpReg, DownReg, Unchanged

    def ExcelDataQuality(inputname, scrappy, repnum):

        ##DataQuality
        print("Producing Data Quality File...")

        subdirectory = "/DataAnalysis/"
        os.makedirs(DayandTime + subdirectory)

        for index, name in enumerate(inputname):

            listheadersSpC = list(str(name) + " SpC R" + str(x + 1) for x in range(0, repnum))
            listheaderslnNSAF = list(str(name) + " lnNSAF R" + str(x + 1) for x in range(0, repnum))
            listheadersSAF = list(str(name) + " SAF R" + str(x + 1) for x in range(0, repnum))
            listheadersloge = list(str(name) + " -log(e) R" + str(x + 1) for x in range(0, repnum))
            combinedheader = ["Protein", "MW", *listheadersSpC, *listheadersloge, *listheadersSAF, *listheaderslnNSAF]
            newpath = os.path.join(DayandTime + subdirectory, str(name) + ".csv")

            with open(newpath, 'w') as myfile:
                wr = csv.writer(myfile, delimiter=',', lineterminator='\n')
                wr.writerow(combinedheader)
                for protein in scrappy[index]:
                    SpC = list(s for s in scrappy[index][protein][0])
                    mw = scrappy[index][protein][1][0]
                    loge = list(n for n in scrappy[index][protein][2])
                    SAF = list(l for l in scrappy[index][protein][3])
                    lnNSAF = list(q for q in scrappy[index][protein][4])
                    combinedrow = [protein, mw, *SpC, *loge, *SAF, *lnNSAF]
                    wr.writerow(combinedrow)

    def GraphValues(filter, name, unfilteredstart, TTest, uniquea, uniqueb):

        print("Producing graphs and figures...")

        subdirectory = "/DataAnalysis/"

        if len(name) == 2:
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

        for y in range(len(name)):
            ig, ax = plt.subplots(1, DesiredReplicateNumber, tight_layout=True, figsize=(16, 4))
            for q in range(DesiredReplicateNumber):
                values = [filter[y][x][2][q] for x in filter[y].keys()]
                N, bins, patches = ax[q].hist(values, bins=50, range=[min(values),max(values)], label=str(name[y]))
                ax[q].title.set_text("Score Values for " + str(name[y]) + " R" + str(q+1))
                fracs = N.astype(float) / N.max()
                norm = colors.Normalize(fracs.min(), fracs.max())
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.spring(norm(thisfrac))
                    thispatch.set_facecolor(color)
            newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " Score Values for Scrappy Data.jpg")
            plt.savefig(newpath)

        for y in range(len(name)):
            ig, ax = plt.subplots(1, DesiredReplicateNumber, tight_layout=True, figsize=(16, 4))
            for q in range(DesiredReplicateNumber):
                values = [filter[y][x][4][q] for x in filter[y].keys()]
                N, bins, patches = ax[q].hist(values, bins=100, range=[min(values),max(values)], label=str(name[y]))
                ax[q].title.set_text("lnNSAF Values for " + str(name[y]) + " R" + str(q+1))
                fracs = N.astype(float) / N.max()
                norm = colors.Normalize(fracs.min(), fracs.max())
                for thisfrac, thispatch in zip(fracs, patches):
                    color = plt.cm.plasma(norm(thisfrac))
                    thispatch.set_facecolor(color)
            newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " lnNSAF Values for Scrappy Data.jpg")
            plt.savefig(newpath)
            plt.close()

    if len(CSVInputName) > 1:
        ExcelDataQuality(CSVInputName, ScrappyIDsDict, DesiredReplicateNumber)
        PrincipleComponent(CSVInputName, ScrappyIDsDict, "Data",0,0, DesiredReplicateNumber)
        TTest, UniqueA, UniqueB = DoTTest(ScrappyIDsDict, Control, Wild, -1, CSVInputName)
        ExcelSampleOutputs(CSVInputName, TTest, UniqueA, UniqueB, -1, 0, Control, Wild, SameSamePValue, DesiredReplicateNumber)
        if len(CSVInputName) == 2:
            PrincipleComponent(CSVInputName, TTest, "TTest", 0, Wild, DesiredReplicateNumber)
            Volcano(CSVInputName, TTest, 0, Wild)
        else:
            PrincipleComponent(CSVInputName, TTest, "TTest", Control, Wild, DesiredReplicateNumber)
            Volcano(CSVInputName, TTest, Control, Wild)
        #HeatMap(CSVInputName, TTest)

    if len(CSVInputName) == 1:
        ExcelDataQuality(CSVInputName, ScrappyIDsDict, DesiredReplicateNumber)
        PrincipleComponent(CSVInputName, ScrappyIDsDict, "Data",0,0, DesiredReplicateNumber)

    GraphValues(ScrappyIDsDict, CSVInputName, AllDataDict, TTest, UniqueA, UniqueB)

    def SameSameAnalysis(replicates, filtered, name):

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
                    combination = [x for x in itertools.combinations(filtered[y][protein][4], 3)]
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
                ax.set_xlabel("Q Value", fontsize=14)
                ax.set_ylabel("% Differentially Regulated IDs", fontsize=14)
                ax.text(float(bisect.bisect_left(chartvalues, 1.0)*0.05) - 0.2,10,str(bisect.bisect_left(chartvalues, 1.0)*0.05) + " ~ 1% PQFDR")
                newpath = os.path.join(DayandTime + subdirectory, str(name[y]) + " Average Q value FDR.jpg")
                plt.savefig(newpath)

                Qvalue1percent.append(bisect.bisect_left(chartvalues, 1.0)*0.05)

            return mean(Qvalue1percent)

    def PvsQGraph(name, UpReg, DownReg, Unchanged, control, wild):

        subdirectory = "/PvQ-Histograms"
        os.makedirs(DayandTime + subdirectory)

        print("Producing P vs Q value venn diagrams")

        if len(name) >= 3:

            for x in range(len(wild)):

                w = int((len(UpReg) / 2) + x)

                plt.clf()
                ##Up
                newpath = os.path.join(DayandTime + subdirectory, "P vs Q Upregulated " + str(name[control]) + " vs " + str(wild[x]) + ".jpg")
                graph = venn2([set(UpReg[x]), set(UpReg[w])], ("P", "Q"))
                graph.get_patch_by_id('10').set_color('blue')
                graph.get_patch_by_id('01').set_color('pink')
                plt.savefig(newpath)

                plt.clf()
                # Down
                newpath = os.path.join(DayandTime + subdirectory, "P vs Q Downregulated " + str(name[control]) + " vs " + str(wild[x]) + ".jpg")
                graph = venn2([set(DownReg[x]), set(DownReg[w])], ("P", "Q"))
                graph.get_patch_by_id('10').set_color('blue')
                graph.get_patch_by_id('01').set_color('pink')
                plt.savefig(newpath)

                plt.clf()
                # Unchanged
                newpath = os.path.join(DayandTime + subdirectory, "P vs Q Unchanged " + str(name[control]) + " vs " + str(wild[x]) + ".jpg")
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
        QValue = SameSameAnalysis(DesiredReplicateNumber, ScrappyIDsDict, CSVInputName)
        print("Q val: " + str(QValue))
        if len(CSVInputName) > 1:
            QTest, UniqueA, UniqueB = DoTTest(ScrappyIDsDict, Control, Wild, QValue, CSVInputName)
            UpReg, DownReg, Unchanged = ExcelSampleOutputs(CSVInputName, QTest, UniqueA, UniqueB, 1, QValue, Control, Wild, SameSamePValue, DesiredReplicateNumber)
            PvsQGraph(CSVInputName, UpReg, DownReg, Unchanged, Control, Wild)

    #End Processes
    def End():
        endtime = time.perf_counter()
        print("Total Processing time: " + str(endtime - starttime))
        shutil.make_archive(os.path.dirname(os.path.realpath(__file__)) + "/Output " + str(CSVInputName), 'zip', DayandTime)
        shutil.rmtree(DayandTime)
    End()

if __name__ == "__main__":
    print("Python Version: " + str(sys.version))
    print("Checking for offline mode")
    if argv != int:
        main("offline")
    else:
        main("online")