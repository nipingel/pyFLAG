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
import time
from multiprocessing import Pool
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
progress bar function 
"""
def progressBar(value, endvalue, beam, bar_length=20):
    beamName = str(beam)
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\rPercent of integrations filled in beam "+ beamName +": [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()  

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
function to sort the input data arrays of shape bank X int X channel into a combined
global data buffer of shape int X full channel range, where full channel range
is either 25*20 = 500 or 20*160 = 3200 nased on mode
"""
def bandpassSort(dataXList, dataYList, xIdList):

    ## determine the number of banks
    numBanks = len(xIdList)
    ## determine size of input arrays
    numInts = dataXList[0].shape[0]
    numChans = dataXList[0].shape[1]

    ## create arrays to hold data
    sortedDataBuff_X = np.zeros([numInts, numChans * 20], dtype = 'float')
    sortedDataBuff_Y = np.zeros([numInts, numChans * 20], dtype = 'float')

    ## loop over array containing the xId and sort based on channel mode
    for bankIdx in range(0, numBanks):
        xID = xIdList[bankIdx]
        ## determine mode; if 25 we are in coarse channel mode,
        ## which requires sorting in a non-contiguous fashion
        if numChans == 25:
            
            ## position in bandpass dictated by xID
            ## starting position in channel space is
            ## xID*5
            bandpassStartChan = xID * 5
            bandpassEndChan = bandpassStartChan + 5

            ## get each chunk of five contigious channels from BANK data,
            ## and place it in the proper spot in full bandpass
            for i in range(0, 5):
                bankStartChan = i * 5
                bankEndChan = bankStartChan + 5
                try:
                    sortedDataBuff_X[:, bandpassStartChan:bandpassEndChan] = dataXList[bankIdx][:, bankStartChan:bankEndChan]
                    sortedDataBuff_Y[:, bandpassStartChan:bandpassEndChan] = dataYList[bankIdx][:, bankStartChan:bankEndChan]
            
                ## in early observations, sometimes an extra integration was recorded in some banks. This pass statement
                ## will skip any extra integrations 
                except IndexError:
                    pass

                ## if a bank stalled, the integrations may not be equal to defined array size. Just skip 
                except ValueError:
                    pass
                ## increment bandpassStartChan/bandpassEndChan by 100 for proper position in full bandpass
                ## recall each chunk of five coarse channels are spaced by 100 in this mode (xId = 0; chans: 0-4, 100-104, ...)
                bandpassStartChan += 100
                bandpassEndChan = bandpassStartChan + 5 
        
        ## if 160 we are in linear fine 
        ## frequency mode. Channels can be
        ## sorted in linear fashion   
        elif numChans == 160:

            ## 160 fine channels are contiguous (e.g., xId = 1; chans:160-319)
            bandpassStartChan = xID * 160
            bandpassEndChan = bandpassStartChan + 160 
            
            ## try-catch to skip over stalled banks
            try:
                sortedDataBuff_X[:, bandpassStartChan:bandpassEndChan] = dataXList[bankIdx][:, :]
                sortedDataBuff_Y[:, bandpassStartChan:bandpassEndChan] = dataYList[bankIdx][:, :]
            except ValueError:
                pass

    return sortedDataBuff_X, sortedDataBuff_Y      

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
function that submits list of bank FITS files to the beamformer module for processing
using multiprocessing's Pool object. Each Bank file is sent to the beamformer object
using one of the available 32 cores on flag-beam.
"""
def multiprocessBankList(bf, bankList, bm):
    p = Pool()
    ## create zipped list to pass in the arguments for each BANK file
    ## beam will always be the same
    processList = list(zip(bankList, bm * len(bankList)))

    """
    start the jobs; result will be a list wherein the
    first element is acess to the XX/YY beamformed spectra (in order of input);
    that is, the first list element is always bank A and last is T. The second element
    gives access to the XX (0-th index)/YY (index 1) spectra. For now, bank label returned
    as index 2. 
    """
    result = p.starmap_async(bf.getSpectralArray, processList)

    ## release worker pool
    p.close()
    p.join()

    ## parse result
    result = np.array(result.get())

    dataXList  = []
    dataYList = []
    dataXArr = np.zeros([len(bankList), result[0][0].shape[0], result[0][0].shape[1]])
    dataYArr = np.zeros([len(bankList), result[0][0].shape[0], result[0][0].shape[1]])
    xIDReturnArray = np.zeros([len(bankList)], dtype='int')
    for idx in range(0, len(bankList)):
        ### ignore if no integrations (stall)
        #if result[idx][0].shape[0] > 0:
            #dataXArr[idx, :, :] = result[idx][0]
            #dataYArr[idx, :, :] = result[idx][1]
        dataXList.append(result[idx][0])
        dataYList.append(result[idx][1])
        xIDReturnArray[idx] = result[idx][2]
    sortedDataBuff_X , sortedDataBuff_Y= bandpassSort(dataXList, dataYList, xIDReturnArray)
    return sortedDataBuff_X, sortedDataBuff_Y


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

    ## split project path to get project name string
    projectPathSplit = projectPath.split('/')
    projectStr = projectPathSplit[-1]
    dataPath = '/mnt/flag/' +  projectStr + '/BF/'

    ##Paths to ancillary FITS files
    goFitsPath = projectPath + '/GO'

    ## test validity of input path
    if not os.path.exists(projectPath):
        print('Incorrect path to project directory; exiting...')
        sys.exit()

    ## get weight FITS files and test validitiy of the path
    wtFileList = glob.glob(weightPath)
    if len(wtFileList) == 0:
        print('Incorrect path to weight FITS files; exiting...')
        sys.exit()

    ## generate object and list of FITS files in case user wishes to process all object/timestamps
    allObjList, allFitsList = generateObjAndFitsLists(goFitsPath)
    
    """
    the existence of optional inputlists of good/bad timestamps (-g/-b flags), objects (-o flag), and beams (-m flag) needs to 
    determined 
    """
    
    ## check if bad time stamps were provided  (-b flag)
    if args.badTimes:
        badScanList = args.badTimes

        ## loop through and append '.fits' to each element for later ID
        badScanList = [s + '.fits' for s in badScanList]
                
        print('\n')
        print('Will remove following timestamps from processing: ')
        print('\n'.join('{}'.format(item) for item in badScanList))

    ## check if explicit list of time stamps were provided (-g flag)
    if args.goodTimes:
        userFitsList = args.goodTimes
        print('\n')
        print('Processing following time stamps: ') 
        print('\n'.join('{}'.format(item) for item in userFitsList))

    ## check if list of objects were provided (-o flag)
    if args.objectList:
        objList = args.objectList
        """
        if the user provides a list of particular objects, first test that the source was observed; if not, 
        program exits. If the user did not provide an object list, process all observed objects. 
        """
        for obj in objList: ## loop 
            if not obj in allObjList:
                print('Requested object, ' + obj + 'was not observed in this session. Exiting...')
                sys.exit()
        print('\n')
        print('Processing following observed objects: ')
        print('\n'.join('{}'.format(item) for item in objList))       

    ## check if list of beams were provided (-m flag)
    if args.beamList:
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
        objList = allObjList
        print('\n')
        print('Processing all observed objects in project...')
    if not args.beamList:
        """
        try-except to define number of beams; if beamList not defined, then 
        get from first weight file header
        """
        try:
            numBeams = int(len(beamList)) 
        except NameError:
            wtFile = fits.open(wtFileList[0])
            numFields = wtFile[1].header['TFIELDS']
            numBeams = int((numFields - 2) / 2)
            ## create list of beam numbers of type str
            beamList = [str(i) for i in range(numBeams)]
            print('\n')
            print('Processing all available beams...')  
    
    ## define beamformingModuleObject
    bf = BeamformingModule(dataPath, weightPath) ## creates beamforming module

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
        if no list was provided, continue processifng all timestamps
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

            """
            Loop over time stamps associated with current observed object
            and process
            """
            fileCnt = 0
            for dataFITSFile in fileList:
                ## inform user of progress
                progressBar(fileCnt, len(fileList), bm, bar_length=20)

                if dataFITSFile[-5:] == '.fits':
                    dataFITSFile = dataFITSFile[:-5]

                ## check if we have bank data for this scan. If not, likely scan aborted early
                testLst = glob.glob(dataPath + dataFITSFile + '*.fits')
                if len(testLst) == 0:
                    print('\n')
                    print('No associated data file with time stamp: ' + dataFITSFile)
                    text_file = open("skippedScans.txt", "w")
                    text_file.write(dataFITSFile)
                    text_file.close()
                    continue

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
                
                sortedDataBuff_X, sortedDataBuff_Y = multiprocessBankList(bf, bankList, bm)
            

                ## add to global data buffers
                globalDataBuff_X_List.append(sortedDataBuff_X)
                globalDataBuff_Y_List.append(sortedDataBuff_Y)
                allBanksList.extend(bankList)
                numBanksList.append(len(bankList))

                ## increment file counter
                fileCnt += 1

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
