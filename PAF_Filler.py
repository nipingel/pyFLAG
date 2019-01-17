# coding: utf-8
"""
1/15/18
This is the main function for the FLAG spectral line filler. It imports the I/O functionality of astropy,
numpy, and matplotlib, gathers the necessary metadata from the GBT ancillary FITS files (e.g. Antenna, GO), and
calls on external modules to perform the beamforming. It then collates these metadata and raw beamformed spectra into a new SDFITS format that is GBTIDL friendly.
User Inputs:
** Note that all list inputs should be space delimited within single quotes
/path/to/project/directory - path to where ancillary FITS files are (e.g. /home/gbtdata/AGBT16B_400)
/path/to/weight/FITS/files - path to weights FITS files; recommended to place '*' wild card in the place of specific bank letter identifier 
restfreq - Rest frequency in Hz (may be phased out once M&C can communicate with IF manager, but now necessary for Doppler corrections)
-b 'List of bad timestamps' - (optional; default not defined) a list of bad timestamps the program should ignore (e.g. '2018_01_15_2017_00:00:00' '2018_01_15_2017_00:00:01') 
-o 'List of objects' - a list of objects to process (e.g. '3C147 NGC6946’); defaults to all objects contained within the observational session.
-g 'List of specific time stamps' - a list of specific timestamps to process (e.g. '2018_01_15_2017_00:00:00' '2018_01_15_2017_00:00:01’); defaults to all time stamps associated with particular observed objects.
-m 'List of beams' -  a list of beams to process (e.g. '1 2 3 4 5 6’); defaults to all beams found in associated weight files.
A single FITS file is produced for each processed beam and is output to which ever directory the call to PAF_Filler is invoked.
Usage:
ipython PAF_Filler.py /path/to/project/ /path/to/weight/FITS/files -b 'List of bad timestamps' -o 'List of objects' -g 'List of specific time stamps' -m 'List of beams'
Eample:
ipython PAF_Filler.py /home/gbtdata/AGBT17B_360_01 /lustre/project/flag/AGBT17B_360_01/BF/weight_files/w_178_02_*.FITS 1420.405e6 -b 2018_08_01:00:00 2018_08_01:52:00 -o NGC891  -g 2018_08_01:05:00 2018_08_01:06:00 2018_08_01:07:00 -m 1 3
__email__ = "nipingel@mix.wvu.edu"
__status__ = "Production"
"""

#imports
from astropy.io import fits
import numpy as np
import sys
import os
import glob
import pickle
from modules.metaDataModule import MetaDataModule
from modules.beamformingModule import BeamformingModule
import matplotlib.pyplot as pyplot

## global variables
numGPU = 10
numBanks = numGPU * 2
pwd = os.getcwd()
pfb = False

## make beam dictionary to map from BYU to WVU convention (e.g. BYU 0 -> WVU 1)
wvuBeamDict = {'0':'1', '1':'2', '2':'6', '3':'0', '4':'3', '5':'5','6':'4'}
byuBeamDict = {'1':'0', '2':'1', '6':'2', '0':'3', '3':'4', '5':'5', '4':'6'}

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
                dataBuff[ints, bandpassStartChan:bandpassEndChan] = bankData[ints, bankStartChan:bankEndChan]
                
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
function to get number of integrations, integration length, and number of channels. Inputs identify the current file
to open. 
"""
def getScanInfo(fitsLst):
    fitsLst.sort() ## sort to be in numerical order before returning
    hdu = fits.open(fitsLst[0]) 
    intLen = np.float(hdu[0].header['REQSTI']) ## integration length
    numInts = np.int(hdu[1].header['NAXIS2']) ## total number of integrations
    form = np.float(hdu[1].header['TFORM3'][:-1]) 
    data = hdu[1].data['DATA']
    numChans = data.shape[1]/2112 ## always 2112 correlations in raw FITS files
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
input, the program processes each requested object in the GBT project by collating all necessary data and ancillary 
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

    ## test for correct number of command line arguments:
    if len(sys.argv) < 3:
        print('Incorrect number of inputs; usage: ->ipython PAF_Filler.py /path/to/project/ /path/to/weight/FITS/files -b \'List of bad timestamps\' \
            -o \'List of objects\' -g \'List of specific time stamps\' -m \'List of beams\'')

        sys.exit()
    
    ## command line inputs
    projectPath = sys.argv[1] ## of the form /home/gbtdata/AGBT16B_400_01
    weightPath = np.str(sys.argv[2]) ## of the form /lustre/projects/flag.old/AGBT16B_400_01/BF/weight_files/*.FITS
    restfreq = np.float(sys.argv[3]) ## in units of Hz
    ## check if an end '/' was given; if so, remove it
    if projectPath[-1] == '/':
        projectPath = projectPath[:-1]

    ## test validity of input path
    if not os.path.exists(projectPath):
        print('Incorrect path to project directory; exiting...')
        sys.exit()

    ## test validitiy of the path to the weight FITS files
    testLst = glob.glob(weightPath)
    if len(testLst) == 0:
        print('Incorrect path to weight FITS files; exiting...')
        sys.exit()
    
    """
    the existance lists of good/bad timestamps (-g/-b flags), objects (-o flag), and beams (-m flag) needs to 
    determinded with if none are provided. 
    """
    ## define list of all flags to compare to later
    totalflagList = ['-b', '-g', '-o', '-m']

    ## define list to append active flags to
    currFlags = []

    ## get index of -b flag
    indB = findIndex('-b', sys.argv)
    ## append to currFlags list if not None type
    if not indB == None:
        currFlags.append('-b')

    ## get index of -g flag
    indG = findIndex('-g', sys.argv)
    ## append to currFlags list if not None type
    if not indG == None:
        currFlags.append('-g')

    ## get index of -o flag 
    indO = findIndex('-o', sys.argv)
    ## append to currFlags list if not None type
    if not indO == None:
        currFlags.append('-o')

    ## get index of -m flag
    indM = findIndex('-m', sys.argv)
    ## append to currFlags list if not None type
    if not indM == None:
        currFlags.append('-m')

    ## place all flags indices in a list
    indList = [indB, indG, indO, indM]

    ## clean list of indices of None type
    indList = [x for x in indList if x is not None]
    indList.sort() ## make numeric order so the numeric order is know 

    ## if indList is not empty, we have some flags to process... 
    if not len(indList) == 0:
        ## loop through list of flags, collect proceeding list elements between current flag and next one. 
        ## deal with flag for bad scan timestamps
        for i, ind in enumerate(indList):
            
            ## process if flag is '-b'
            if ind == indB:
                """
                the bad scans will be all elements of the list between this index and the next in indList.                
                if there is only one element of indList (only one flag was set), or we are on the last element of indList, 
                just take remaining elments from sys.argv list
                """
                if i == len(indList) - 1:
                    badScanList = sys.argv[ind + 1:]
                else:
                    badScanList = sys.argv[ind + 1:indList[i + 1 ]]

                ## loop through and append '.fits' to each element for later ID
                badScanList = [s + '.fits' for s in badScanList]
                
                ## inform user
                if len(badScanList) == 0:
                    print('\n')
                    print('No bad time stamp list after \'-b\' flag; exiting...')
                    sys.exit()
                else:
                    print('\n')
                    print('Will remove following timestamps from processing: ')
                    print('\n'.join('{}'.format(item) for item in badScanList))

            ## process if flag is '-g'
            elif ind == indG:
                if i == len(indList) - 1:
                    userFitsList = sys.argv[ind + 1:]
                else:
                    userFitsList = sys.argv[ind + 1:indList[i + 1 ]]

                
                ## inform user
                if len(userFitsList) == 0:
                    print('\n')
                    print('No time stamp list after \'-g\' flag; exiting...')
                    sys.exit()
                else:
                    print('\n')
                    print('Processing following time stamps: ') 
                    print('\n'.join('{}'.format(item) for item in userFitsList))

            ## process if flag '-o'
            elif ind == indO:
                if i == len(indList) - 1:
                    obsObjectList = sys.argv[ind + 1:]
                else:
                    obsObjectList = sys.argv[ind + 1:indList[i + 1 ]]
                
                ## inform user
                if len(obsObjectList) == 0:
                    print('\n')
                    print('No observed objects list after \'-o\' flag; exiting...')
                    sys.exit()
                else:
                    print('\n')
                    print('Processing following observed objects: ')
                    print('\n'.join('{}'.format(item) for item in obsObjectList))

            elif ind == indM:
                if i == len(indList) - 1:
                    beamList = sys.argv[ind + 1:]
                else:
                    beamList = sys.argv[ind + 1:indList[i + 1 ]]
                ## loop through beamList to translate from WVU convention to BYU for consistent processing
                #byuBeamList = []
                #for bm in beamList:
                #   byuBeamList.append(byuBeamDict[bm])

                ## inform user
                if len(beamList) == 0:
                    print('\n')
                    print('No beam list after \'-m\' flag; exiting...')
                    sys.exit()
                else:
                    print('\n')
                    print('Processing following beams: ')
                    print('\n'.join('{}'.format(bm) for bm in beamList))
            ## if none of these flags were found, then inform user that an unrecognized flag was given and exit
            else:
                print("Unrecognized flag was input. Accepted flags are: -b, -g, -o, and -m. Exiting...")
                sys.exit()
    ## one more if block to inform user if any default settings active
    if not '-b' in currFlags:
        print('\n')
        print('No bad time stamps indicated...')
    if not '-g' in currFlags:
        print('\n')
        print('Processing all available scans...')
    if not '-o' in currFlags:
        print('\n')
        print('Processing all observed objects in project...')
    if not '-m' in currFlags:
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
        #numBeams = len(byuBeamList)   
        numBeams = len(beamList) 
    except NameError:
        numFields = wtFile[1].header['TFIELDS']
        numBeams = (numFields - 2) / 2 
        #byuBeamList = range(0, numBeams)
        beamList = range(0, numBeams)
        #for el in range(0, len(byuBeamList)):
        #   byuBeamList[el] = np.str(el)
        for el in range(0, len(beamList)):
           beamList[el] = np.str(el)
 
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
	#fileList = fileList[3:45] ## AGBT16B_400_12, HI Grid
        #fileList = fileList[31:85] ## AGBT16B_400_12, 3C295 Grid
        #fileList = fileList[0:41] ## AGBT16B_400_13, HI Grid
        #fileList = fileList[1] + fileList[10] ## AGBT16B_400_14 3C147 7Pt Cal Scan
        #fileList= fileList[0:96] ## AGBT17B_455_01 -> first 7Pt-Cal Scan
	fileList = fileList[1:]
        #fileList = fileList[105:-5]
	#fileList = fileList[3:117]
	#fileList = fileList[117:]
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
        #for bm in byuBeamList:
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
            if isinstance(fileList, unicode):
                 fileList = [fileList]
            """
            Loop over time stamps associated with current observed object 
            """
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
                    print('Skipping to next scan...')
                    text_file = open("skippedScans.txt", "w")
                    text_file.write(dataFITSFile)
                    text_file.close()
                    continue 
                else:
                    ## get essential info about current timestamp
                    numInts, intLen, numSpecChans, bankList = getScanInfo(testLst)

                    ## check whether we are in PFB (fine channelization) mode. If so, set pfb to True
                    if numSpecChans == 160:
                        pfb = True
                    else:
                        pfb = False
                    """
                    initialize bank data buffers; columrns are total number of spectral channels and rows are 
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
            build metadata; inputs are path to ancillary FITS files, path to raw BANK FITS files, names of the processed BANK files for this beam
            the global data buffers holding the beamformed spectra, the beam that was processed, and the boolean denoting which spectral mode we're in.          
            """
            md = MetaDataModule(projectPath, dataPath, weightPath, fileList, restfreq, allBanksList, numBanksList, globalDataBuff_X_List, globalDataBuff_Y_List, bm, pfb)
            thduList = md.constuctBinTableHeader()
            dataFITSFile = dataFITSFile[:-6]
	    ## save the file
            thduList.writeto(pwd+'/' + projectStr + '_' + source + '_Beam'+str(bm) + '.fits') ## finally, saves the file.
    
       
    
#run PAF_Filler.py (main) function
if __name__=="__main__":
    main()
