#### coding: utf-8
"""
10/31/19
This is the contains the main function for the FLAG spectral line filler. It gathers the necessary metadata from 
the GBT ancillary FITS files (e.g. Antenna, GO), and calls on external modules to perform the beamforming. The beamforming
module is called in a multi-process fashion where each core handles a BANK file. It then collates these metadata and raw 
beamformed spectra into a new SDFITS format that is GBTIDL friendly.

User Inputs:

/path/to/project/directory - <required> path to where ancillary FITS files are (e.g. /home/gbtdata/AGBT16B_400)
/path/to/weight/FITS/files - <required> path to weights FITS files; recommended to place '*' wild card in the place of specific bank letter identifier 
restfreq - <required> Rest frequency in Hz (may be phased out once M&C can communicate with IF manager, but now necessary for Doppler corrections)
centralFreq - <required> central frequency in Hz (may be phased out once M&C can communicate with IF manager, but now necessary when LO is taken out of scan coordinator)
-b --badTimes -<optional> a list of bad timestamps the program should ignore (e.g. '2018_01_15_2017_00:00:00' '2018_01_15_2017_00:00:01'); no default
-o --objectList - <optional> a list of objects to process (e.g. '3C147 NGC6946’); defaults to all objects contained within the observational session.
-g  --goodTimes - <optional> a list of specific timestamps to process (e.g. '2018_01_15_2017_00:00:00' '2018_01_15_2017_00:00:01’); defaults to all time stamps associated with particular observed objects.
-m --beamList - <optional> a list of beams to process (e.g. '1 2 3 4 5 6’); defaults to all beams found in associated weight files.
A single FITS file is produced for each processed beam and is output to which ever directory the call to PAF_Filler is invoked.
Usage:
ipython --c="run PAF_Filler.py /path/to/project/ /path/to/weight/FITS/files restfreq [Hz] centralFreq [Hz] -b 'List of bad timestamps' -o 'List of objects' -g 'List of specific time stamps' -m 'List of beams'""
Eample:
ipython --c="run PAF_Filler.py /home/gbtdata/AGBT17B_360_01 /lustre/project/flag/AGBT17B_360_01/BF/weight_files/w_178_02_*.FITS 1420.405e6 1450e6 -b 2018_08_01:00:00 2018_08_01:52:00 -o NGC891  -g 2018_08_01:05:00 2018_08_01:06:00 2018_08_01:07:00 -m 1 3"
__email__ = "Nickolas.Pingel@anu.edu.au"
__status__ = "Production"
"""

#imports
from astropy.io import fits
import numpy as np
import sys
import os
import glob
import pickle
import argparse
from modules.metaDataModule import MetaDataModule
from modules.beamformingModule import BeamformingModule
import matplotlib.pyplot as pyplot

## global variables
numGPU = 10
numBanks = numGPU * 2
pwd = os.getcwd()
pfb = False

## beam dictionaries to map from between BYU to WVU beam name conventions (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {'0':'1', '1':'2', '2':'6', '3':'0', '4':'3', '5':'5','6':'4'}
byuBeamDict = {'1':'0', '2':'1', '6':'2', '0':'3', '3':'4', '5':'5', '4':'6'}

"""
Takes input lists and sorts them into a 'list-of-lists' that contains the 
original elements sorted into different observing modes. Each element in the
returned list-of-lists is the portion of the input lists associated with a 
particular observing mode.
"""
def sortLists(dataXList, dataYList, allFilesList, allBanksList, allNumBanksList):
    """
    define lists to store the sorted input lists if needed
    """
    allSortedDataX = []
    allSortedDataY = []
    allSortedFiles = []
    allSortedBanks = []
    allSortedNumBanks = []

    """
    loop through the global data buffers to determine the mode based on spectral channels 
    """
    modeSwIdx = 0
    for idx in range(0, len(dataXList)):
        dataArr = dataXList[idx]
        """
        If first iteration, set the prevNumChans to value of first element. 
        Recall that the numpy arrays are of shape [ints, specChans]
        """
        if idx == 0:
            prevNumChans = dataArr.shape[1]

        """
        Test whether current channel length is equal to previous; 
        if not, the mode has changed
        """
        currNumChans = dataArr.shape[1]
        if prevNumChans != currNumChans:
            """
            Append to sorted lists data from the last mode switch idx to
            current idx (python indices to one less)
            """
            allSortedDataX.append(dataXList[modeSwIdx:idx])
            allSortedDataY.append(dataYList[modeSwIdx:idx])
            allSortedFiles.append(allFilesList[modeSwIdx:idx])
            allSortedNumBanks.append(allNumBanksList[modeSwIdx:idx])
            
            """
            The total number of list elements to add to allSortedBanks
            list is equal the sum of allNumBanksList elements
            between 0 and current (if first mode change), or 
            --- if subsequent iteration --- the sum of allNumBanksList 
            elements from 0 to the index of last mode change to the sum 
            of allNumBanksList elements to current index. 
            """
            if modeSwIdx == 0:
                allSortedBanks.append(allBanksList[0:int(np.sum(allNumBanksList[0:idx]))])
            else:
                allSortedBanks.append(allBanksList[int(np.sum(allNumBanksList[0:modeSwIdx])):int(np.sum(allNumBanksList[0:modeSwIdx])) + int(np.sum(allNumBanksList[modeSwIdx:idx]))])
            modeSwIdx = idx

        """
        upon the last iteration, append everything to sorted lists from
        last mode switch idx (could still be zero) to current idx + 1 to 
        ensure last element is not left behind
        """
        if idx == len(dataXList) -1:
            allSortedDataX.append(dataXList[modeSwIdx:idx + 1])
            allSortedDataY.append(dataYList[modeSwIdx:idx + 1])
            allSortedFiles.append(allFilesList[modeSwIdx:idx + 1])
            allSortedNumBanks.append(allNumBanksList[modeSwIdx:idx + 1])
            if modeSwIdx == 0 and idx == 0:
                allSortedBanks.append(allBanksList[0:int(np.sum(allNumBanksList[0:idx + 1]))])
            else:
                allSortedBanks.append(allBanksList[int(np.sum(allNumBanksList[0:modeSwIdx])):int(np.sum(allNumBanksList[0:modeSwIdx])) + int(np.sum(allNumBanksList[modeSwIdx:idx + 1]))])
        else:
            prevNumChans = currNumChans
    
    ## return sorted lists
    return allSortedDataX, allSortedDataY, allSortedFiles, allSortedBanks, allSortedNumBanks

"""
function to map individual BANK data (bankData)
to the full bandpass data buffer (dataBuff) based on XID. 
"""
def bandpassSort(xID, dataBuff, bankData):
    """
    Determine the correlation mode by 
    the number of channels stored in
    data buff as the second element.
    Note that the number of integrations dataBuff's first element
    """
    numInts = bankData.shape[0]
    numChans = bankData.shape[1]
    for ints in range(0, numInts):
        if numChans == 25:
            
            ## coarse channel mode
            ## position in bandpass dictated by xID
            bandpassStartChan = xID * 5
            bandpassEndChan = bandpassStartChan + 5
            
            ## get each chunk of five contigious channels from BANK data,
            ## and place it in the proper spot in full bandpass
            for i in range(0, 5):
                bankStartChan = i * 5
                bankEndChan = bankStartChan + 5
                try:
                    dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, bankStartChan:bankEndChan]
                
                ## in early observations, sometimes an extra integration was recorded in some banks. This pass statement
                ## will skip any extra integrations  
                except IndexError:
                    pass
                ## increment bandpassStartChan/bandpassEndChan by 100 for proper position in full bandpass
                bandpassStartChan += 100
                bandpassEndChan = bandpassStartChan + 5
        elif numChans == 160:
            bandpassStartChan = xID * 160
            bandpassEndChan = bandpassStartChan + 160

            ## Added to avoid error where BANK data length is zero from a stall
            if len(bankData[ints,:]) == 0:
                dataBuff[ints, bandpassStartChan:bandpassEndChan] = np.zeros([160], dtype='float32')
            elif ints >= len(dataBuff[:, 0]):
                continue
            else:    
                dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, :]
    return dataBuff

"""
Function that generates and returns list of observed objects and the associated timestamp FITS files. 
The GO FITS files are searched and sorted to generate the required lists, which is why the path to these 
files is an input to this function. 
"""
def generateObjAndFitsLists(goFitsPath):        

    goFits = glob.glob(goFitsPath+'/*.fits') ## read in GO FITS files to
    goFits.sort() ## sort to get correct time stamp order
    genObjList = [] ## define empty list to be returned (after casting to list for indexing)
    """
    define empty list. Once we add new object, an empty list will be appended and
    associated FITS files added. For example, if the second object added is '3C147', 
    genFitsList[1] will be a list of all associated timestamp FITS files. 
    """
    genFitsList = [] 
    itr = 0 ## counter variable to check if first iteration
    ## iterate over all GO FITS files to check the observed object
    for goName in goFits:
       goHDU = fits.open(goName)
       if itr == 0: ## upon first iteration, append an empty list to the generated list of FITS files
           genObjList.append(goHDU[0].header['OBJECT'])
           genFitsList.append([])
           genFitsList[0].append(goName[-24:])
           itr+=1
       else: ## now processing for remaining 
           obj = goHDU[0].header['OBJECT'] 
           if not any(x == obj for x in genObjList): ## test if object already in list
               genObjList.append(obj) 
               genFitsList.append([]) ## append empty list
               genFitsList[-1].append(goName[-24:])
           else:
            ind = genObjList.index(obj)
            genFitsList[ind].append(goName[-24:])
    return genObjList, genFitsList

"""
function to get number of integrations, integration length, and number of channels and BANK file list. Input identifies the current file
to open. 
"""
def getScanInfo(fitsLst):
    fitsLst.sort() ## sort to be in numerical order before returning
    hdu = fits.open(fitsLst[0]) 
    intLen = np.float(hdu[0].header['REQSTI']) ## integration length
    numInts = np.int(hdu[1].header['NAXIS2']) ## total number of integrations
    form = np.float(hdu[1].header['TFORM3'][:-1]) 
    data = hdu[1].data['DATA']
    numChans = np.int(data.shape[1]/2112) ## always 2112 correlations in raw FITS files
    return numInts, intLen, numChans, fitsLst

"""
small function to return the index of the indicated string in the provided list. This function is used primarly to parse
user input
"""
def findIndex(elemStr, inputList):
    try:
        index = inputList.index(elemStr)
        return index
    except ValueError:
        return None


"""
This is the 'main' function that drives the creation of an SDFITS file for the selected beams. After handling the user 
input, the function processes each requested object in the GBT project by collating all necessary data and ancillary 
FITS files. The beamformer object (bf) is called to perform the beamforming. Several processing functions within this
main script are called to sort the returned beamform spectrum into its final data container. Once the data are correctly
formatted, the metadata object (md) is called to correctly format and sort the metadata from the GO and Antennna ancillary
FITS files.
""" 
def main():

    ## EXCEPTION HANDLING FOR USER INPUT
    """
    since we are writing out FITS files in the users's current directory, we need to check 
    the correct write permissions are in place. Try to open a file and catch an exception
    """
    try:
        tempFileName = pwd + 'testFile'
        f = open(tempFileName, 'w')
        f.close()
        os.remove(tempFileName)
    except PermissionError:
        print('Please move to a directory where you have adaquete write permissions; exiting...')
        sys.exit()

    ## define arguement parser
    parser = argparse.ArgumentParser()

    ## add positional and required arguments
    parser.add_argument("projectPath", help="<required> path to ancillary FITS files (e.g. /home/gbtdata/AGBT16B_400)")
    parser.add_argument("weightPath", help="<required> path to weights FITS files")
    parser.add_argument("restFreq", help="<required> rest frequency in Hz (may be phased out once M&C can communicate with IF manager, but now necessary for Doppler corrections", type = float)
    parser.add_argument("centralFreq", help="<required> central frequency in Hz (may be phased out once M&C can communicate with IF manager, but now necessary when LO is taken out of scan coordinator)", type = float)

    ## add optional short options
    parser.add_argument("-b", "--badTimes", help="<optional> a list of bad timestamps the program should ignore (e.g. -b 2018_01_15_2017_00:00:00 2018_01_15_2017_00:00:01); no default", nargs = "+")
    parser.add_argument("-o", "--objectList", help="<optional> a list of objects to process (e.g. -o 3C147 NGC6946); defaults to all objects contained within the observational session.", nargs = "+")
    parser.add_argument("-g", "--goodTimes", help="<optional> a list of specific timestamps to process (e.g. -g 2018_01_15_2017_00:00:00 2018_01_15_2017_00:00:01); defaults to all time stamps associated with particular observed objects", nargs = "+")
    parser.add_argument("-m", "--beamList", help="<optional> a list of beams to process (e.g. -m 1 2 3 4 5 6); defaults to all beams found in associated weight files.", nargs = "+")

    ## finally, parse arguments
    args, unknown = parser.parse_known_args()

    ## assign required command line inputs
    projectPath = args.projectPath ## of the form /home/gbtdata/AGBT16B_400_01
    weightPath = args.weightPath + '/*.FITS' ## of the form /lustre/projects/flag.old/AGBT16B_400_01/BF/weight_files/*.FITS
    restfreq = args.restFreq ## in units of Hz
    centralfreq = args.centralFreq ## in units of Hz
   
    ## check if an end '/' was given; if so, remove it
    if projectPath[-1] == '/':
        projectPath = projectPath[:-1]

    ## test validity of input path
    #if not os.path.exists(projectPath):
    #    print('Incorrect path to project directory; exiting...')
    #    sys.exit()

    ## test validitiy of the path to the weight FITS files
    #testLst = glob.glob(weightPath)
    #if len(testLst) == 0:
    #    print('Incorrect path to weight FITS files; exiting...')
    #    del testLst
    #    sys.exit()
    
    """
    the existence of optional inputlists of good/bad timestamps (-g/-b flags), objects (-o flag), and beams (-m flag) needs to 
    determined 
    """
    
    ## check if bad time stamps were provided  (-b flag)
    if args.badTimes not None:
        badScanList = args.badTimes

        ## loop through and append '.fits' to each element for later ID
        badScanList = [s + '.fits' for s in badScanList]
                
        print('\n')
        print('Will remove following timestamps from processing: ')
        print('\n'.join('{}'.format(item) for item in badScanList))

    ## check if explicit list of time stamps were provided (-g flag)
    elif args.goodTimes not None:
        userFitsList = args.goodTimes
        print('\n')
        print('Processing following time stamps: ') 
        print('\n'.join('{}'.format(item) for item in userFitsList))

    ## check if list of objects were provided (-o flag
    elif args.objectList not None:
        obsObjectList = args.objectList
        print('\n')
        print('Processing following observed objects: ')
        print('\n'.join('{}'.format(item) for item in obsObjectList))       

    ## check if list of beams were provided (-m flag)
    elif args.beamList not None:
        beamList = args.beamList
        print('\n')
        print('Processing following beams: ')
        print('\n'.join('{}'.format(bm) for bm in beamList))

    ## one more if block to inform user if any default settings are active
    if not args.badTimes:
        print('\n')
        print('No bad time stamps indicated...')
    if not args.goodTimes:
        print('\n')
        print('Processing all available scans...')
    if not args.objectList:
        print('\n')
        print('Processing all observed objects in project...')
    if not args.beamList:
        print('\n')
        print('Processing all available beams...')

    ## split project path to get project name string
    projectPathSplit = projectPath.split('/')
    projectStr = projectPathSplit[-1]
    dataPath = '/lustre/flag/' +  projectStr + '/BF/'

    ##Paths to ancillary FITS files
    goFitsPath = projectPath + '/GO'  
    
    ## define beamformingModuleObject
    bf = BeamformingModule(dataPath, weightPath) ## creates beamforming module
    
    ## open weight file to get number of beams
    wtFileList = glob.glob(weightPath)
    wtFile = fits.open(wtFileList[0])

    """
    try-except to define number of beams; if beamList not defined, then 
    get from first weight file header
    """
    try:
        numBeams = int(len(beamList)) 
    except NameError:
        numFields = wtFile[1].header['TFIELDS']
        numBeams = int((numFields - 2) / 2)
        ## create list of beam numbers of type str
        beamList = [str(i) for i in range(numBeams)]
 
    ## generate object and list of FITS files in case user wishes to process all object/timestamps
    allObjList, allFitsList = generateObjAndFitsLists(goFitsPath)

    """
    if the user has already provided a list of particular objects, first test that the source was observed; if not, 
    program exits. If the user did not provide an object list, process all observed objects. 
    """
    try:
        objList = obsObjectList
        for obj in obsObjectList: ## loop 
            if not obj in allObjList:
                print('Requested object, ' + obj + 'was not observed in this session. Exiting...')
                sys.exit()
    ## if objObjectList is not defined, we are processing all objects
    except NameError: 
        objList = allObjList
        pass

    ## inform user of the inputs and objects to process
    print('\n')
    print('Project path: ' + projectPath)
    print('Data path: ' + dataPath)
    print('Weight file path: ' + weightPath)
    print('Observed objects to process: ')
    print('\n'.join('{}'.format(item) for item in objList))
    

    """
    This is the first of several processing loops. Here we being by looping over the list of objects...
    """
    for objs in objList:
        objInd = allObjList.index(objs)
        source = allObjList[objInd] ## get source
        print('\n')
        print('Processing the source:')
        print(source)
        
        """
        if user has already provided a list of timestamps for processing, test that the timestamps are valid;
        if no list was provided, continue processing all timestamps
        """
        try:
            fileList = userFitsList
            for timeStamp in userFitsList: ## loop
                if not timeStamp + '.fits' in allFitsList[objInd] :
                    print('Requested timestamps, ' + timeStamp + ' is not available for this object.')
                    sys.exit()
        except NameError:
            fileList = allFitsList[objInd]
            pass
	## remove bad scans from file list
        if 'badScanList' in locals():
            for s in badScanList:
                fileList.remove(s)
                try:
                    fileList.remove(s)
                except ValueError:
                    if s == badScanList[-1]:
                        pass
                    else:
                        continue
        print('\n')
        print('Processing scans: ')
        ## loop through and append '.fits' to each element for later ID
        tmpVar = fileList[0]
        if tmpVar[-5:] == '.fits':
             fileList = [s[:-5] for s in fileList] 
        print('\n'.join('{}'.format(item) for item in fileList))
        
        """
        Loop over beams to read in files for object of interest and construct a single SINGLE DISH binary FITS table
        """
        for bm in beamList:

            allBanksList = [] ## list of paths to all good BANK FITS files associated with object
            numBanksList = [] ## number of good BANKS associated with a scan

            ## global lists to hold data. Each element of the list contains a numpy array in the shape of (ints, freq chans)
            globalDataBuff_X_List = []
            globalDataBuff_Y_List = []
                
            """
            to avoid the for loop simply looping over a string (when len of fileList = 1), set if statement
            such that if type of fileList is equal to string, cast into a list
            """
            if isinstance(fileList, str):
                 fileList = [fileList]

            """
            Loop over time stamps associated with current observed object to make sure
            data exists. If not, then append these time stamps to a list for removal
            """
            removeList = []
            for dataFITSFile in fileList:
                """ 
                check if '.fits' is at the end of elements in fileList
                if so, remove it
                """
                if dataFITSFile[-5:] == '.fits':
                    dataFITSFile = dataFITSFile[:-5]
                ## check if we have bank data for this scan
                testLst = glob.glob(dataPath + dataFITSFile + '*.fits')
                if len(testLst) == 0:
                    print('\n')
                    print('No associated data file with time stamp: ' + dataFITSFile)
                    text_file = open("skippedScans.txt", "w")
                    text_file.write(dataFITSFile)
                    text_file.close()
                    removeList.append(dataFITSFile)
            ## remove any scans with no data (probably aborted early)
            if len(removeList) > 0:
                for replaceElem in removeList:
                    fileList.remove(replaceElem)

            """
            Loop over time stamps associated with current observed object
            and process
            """
            for dataFITSFile in fileList:
                if dataFITSFile[-5:] == '.fits':
                    dataFITSFile = dataFITSFile[:-5]
                ## check if we have bank data for this scan
                testLst = glob.glob(dataPath + dataFITSFile + '*.fits')
                ## get essential info about current timestamp
                numInts, intLen, numSpecChans, bankList = getScanInfo(testLst)
                
                """
                initialize bank data buffers; columns are total number of spectral channels and rows are 
                individual integrations
                """
                dataBuff_X = np.zeros([numInts, numSpecChans * numBanks])
                dataBuff_Y = np.zeros([numInts, numSpecChans * numBanks])

                ## initialize weight data buffers in case there is analysis to be done
                xWeightBuff =  np.zeros([numSpecChans * numBanks], dtype = 'complex64')
                yWeightBuff =  np.zeros([numSpecChans * numBanks], dtype = 'complex64')

                bankCnt = 0 ## set bank counter to keep track of number of bank FITS files processed
                """
                Loop through list of BANK FITS files to process 1/20th of the bandpass at a time. Once information 
                about the scan is retrieved, the 'getSpectralArray' in the bf module is called and a beamformed 
                spectrum is returned (rows: inegrations; columns: frequency channels). The weights used in the processing
                are also returned. The function bandpassSort is then called to place this 1/20th chunk of bandpass at the 
                correct position in the global data buffers (that will be written in the output FITS file). 
                """
                for fileName in bankList:
                    print('\n')                
                    print('Beamforming correlations in: '+fileName[-25:]+', Beam: ' + bm) 
                    ## bank name is ALWAYS sixth-to-last character in string
                    bank = fileName[-6] 
                    corrHDU = fits.open(fileName) ## grab xid from dictionary for sorting

                    nRows = corrHDU[1].header['NAXIS2']
                    data = corrHDU[1].data['DATA']
                    xID = np.int(corrHDU[0].header['XID'])

                    if nRows != 0: ## to avoid an empty FITS file
                        """
                        Do the beamforming; returns processed BANK data (cov matrices to a beam-formed bandpass) in both
                        XX/YY Pols; returned data are in form rows: integrations, columns: frequency channels
                        """
                        xData, yData = bf.getSpectralArray(fileName, data, bm, xID)

                        """
                        Sort the 1/20th chunk based on xid number for each integration. The dataBuff_X/Y variable is constantly
                        being updated as the individual BANKS are processed
                        """       
                        dataBuff_X = bandpassSort(xID, dataBuff_X, xData)
                        dataBuff_Y = bandpassSort(xID, dataBuff_Y, yData)
                        allBanksList.append(fileName) ## append good BANK file to master list that is passed to metadataModule
                        bankCnt += 1 ## increment bank counter

                globalDataBuff_X_List.append(dataBuff_X) ## append to global buffer lists
                globalDataBuff_Y_List.append(dataBuff_Y)
                numBanksList.append(bankCnt) ## append bankCnt to list 

            print('\n') ## make terminal output look cleaner

            """
            Because the same object could have been observed in different modes (i.e., CALCORR or PFBCORR), and  
            we wish to have a single SDFITS file per beam, we need to check the global data buffers to change in 
            stored numpy array shape. Send the fileList, global data buffers, allBanksList, and numBanksList to 
            a sorting function. The returned lists contain the data and other necessary inputs to the metadata module
            sorted by operating mode
            contain the data from
            """
            sortedDataBufferList_X, sortedDataBufferList_Y, sortedFileList, sortedBanksList, sortedNumBanksList = sortLists(globalDataBuff_X_List, globalDataBuff_Y_List, fileList, allBanksList, numBanksList)
            """
            Loop over the elements of the sorted lists. sortedDataBuffer_X has a single element, there were no mode changes throughout the observation. 
            In this case, create primary header and write out file. If sortedDataBuffer_X has more than a single element, it means there were multiple mode changes. 
            In this case, create primary header on first iteration and attach the returned binary FITS tables before writing out. 
            """
            for metaIdx in range(0, len(sortedDataBufferList_X)):
                """
                build metadata; inputs are path to ancillary FITS files, path to raw BANK FITS files, names of the processed BANK files for this beam
                the global data buffers holding the beamformed spectra, the beam that was processed, and the boolean denoting which spectral mode we're in.          
                """

                globalDataBuff_X = sortedDataBufferList_X[metaIdx]
                globalDataBuff_Y = sortedDataBufferList_Y[metaIdx]
                sortedFiles = sortedFileList[metaIdx]
                sortedBanks = sortedBanksList[metaIdx]
                sortedNumBanks = sortedNumBanksList[metaIdx]
                
                ## check whether we are in PFB (fine channelization) mode. If so, set pfb to True
                specChans = globalDataBuff_X[0].shape[1]
                if specChans == 3200:
                    pfb = True
                else:
                    pfb = False
                """
                Check if first iteration. If so, set firstIter to True to inform metadata instance
                to create a primary HDU
                """
                if metaIdx == 0:
                    firstIter = True
                else:
                    firstIter = False
                ## instantiate metadata object    
                md = MetaDataModule(projectPath, dataPath, weightPath, sortedFiles, restfreq, centralfreq, sortedBanks, sortedNumBanks, globalDataBuff_X, globalDataBuff_Y, bm, pfb, firstIter)

                if firstIter == True:                
                    ## construct the binary table & header
                    thduList = md.constuctBinTableHeader()
                else:
                    nextThduList = md.constuctBinTableHeader()
                    thduList.append(nextThduList)

	        ## FINALLY, save the file
            thduList.writeto(pwd+'/' + projectStr + '_' + source + '_Beam'+str(bm) + '.fits')
    
#run PAF_Filler.py (main) function
if __name__=="__main__":
    main()
