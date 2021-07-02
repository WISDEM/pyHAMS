import os
import os.path as osp
import platform
import numpy as np

if not platform.system() == 'Windows':
    import pyhams.libhams as hams
else:
    import subprocess as sub

def nemohmesh_to_pnl(nemohMeshPath, writeDir=None):
    '''
    convert mesh from .nemoh format to HAMS .pnl format

    Parameters
    ----------
    nemohMeshPath: str
        path to the nemoh mesh file.

    Returns
    -------
    None

    Raises
    ------
    None

    Notes
    -----
    Nemoh mesh must have:
      - single line header
      - line beginning with '0' separating points/panels
      - line beginning with '0' at the end of the file
      - no duplicate points (nodes)
    '''

    nemohMeshPath = osp.normpath(nemohMeshPath)
    nemohDirName, nemohFileName = osp.split(nemohMeshPath)
    if writeDir is None:
        writeDir = nemohDirName
    if osp.isdir(writeDir) is not True:
        os.makedirs(writeDir)

    # N.B. Nemoh input files can have slightly different headers (panels and
    # points on top line...or just '2 0' or '2 1' on top line).
    iFile = open(nemohMeshPath, 'r')
    lines = iFile.readlines()
    header = lines[0].split()
    if header[0] == '2':
        ySym = int(header[1])
    else:
        ySym = 0
    ixOnes = []
    ixZeros = []
    for ix, line in enumerate(lines):
        if line.split()[0] == '0':
            ixZeros.append(ix)
        if line.split()[0] == '1':
            ixOnes.append(ix)

    numHeaders = ixOnes[0]
    numVertices = ixZeros[0] - numHeaders
    numPanels = ixZeros[1] - 1 - numVertices - numHeaders

    oFilePath = osp.join(writeDir, 'HullMesh.pnl')

    oFile = open(oFilePath, 'w')
    oFile.write('    --------------Hull Mesh File---------------\n\n')
    oFile.write('    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n')
    oFile.write(f'         {numPanels}         {numVertices}         0         {ySym}\n\n')
    oFile.write('    #Start Definition of Node Coordinates     ! node_number   x   y   z\n')

    for ix, line in enumerate(lines[numHeaders:]):
        if line.split()[0] == '0':
            oFile.write('   #End Definition of Node Coordinates\n\n')
            oFile.write('   #Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID)\n')
            break
        oFile.write(f'{line.split()[0]:>5}{line.split()[1]:>18}{line.split()[2]:>18}{line.split()[3]:>18}\n')
    for ix, line in enumerate(lines[(numVertices+2):]):
        if line.split()[0] == '0':
            oFile.write('   #End Definition of Node Relations\n\n')
            oFile.write('    --------------End Hull Mesh File---------------\n')
            break
        if line.split()[0] == line.split()[3]:
            numNodes = 3
            oFile.write(f'{(ix+1):>5}{numNodes:>5}{line.split()[0]:>10}{line.split()[1]:>10}{line.split()[2]:>10}\n')
        else:
            numNodes = 4
            oFile.write(f'{(ix+1):>5}{numNodes:>5}{line.split()[0]:>10}{line.split()[1]:>10}{line.split()[2]:>10}{line.split()[3]:>10}\n')
    oFile.close()


def create_hams_dirs(baseDir=None):
    '''
    create necessary HAMS directories in baseDir

    Parameters
    ----------
    baseDir: str
        The top directory in which to create HAMS Input and Output directories

    Returns
    -------
    None

    Raises
    ------

    Notes
    -----

    '''

    if baseDir is None:
        baseDir = os.getcwd()
    else:
        baseDir = osp.normpath(baseDir)
        if osp.isdir(baseDir) is not True:
            os.makedirs(baseDir)

    inputDir = osp.join(baseDir, 'Input')
    outputDirHams = osp.join(baseDir, 'Output/Hams_format')
    outputDirHydrostar = osp.join(baseDir, 'Output/Hydrostar_format')
    outputDirWamit = osp.join(baseDir, 'Output/Wamit_format')

    if osp.isdir(inputDir) is not True:
        os.mkdir(inputDir)
    if osp.isdir(outputDirHams) is not True:
        os.makedirs(outputDirHams)
    if osp.isdir(outputDirHydrostar) is not True:
        os.makedirs(outputDirHydrostar)
    if osp.isdir(outputDirWamit) is not True:
        os.makedirs(outputDirWamit)

def write_hydrostatic_file(projectDir=None, cog=np.zeros(3), mass=np.zeros((6,6)),
                           dampingLin=np.zeros((6,6)), dampingQuad=np.zeros((6,6)),
                           kHydro=np.zeros((6,6)), kExt=np.zeros((6,6))):
    '''
    Writes Hydrostatic.in for HAMS

    Parameters
    ----------
    projectDir: str
        main HAMS project directory - Hydrostatic.in will be written in ./Input/
    cog: array
        3x1 array - body's CoG
    mass: array
        6x6 array - body's mass matrix
    dampingLin: array
        6x6 array - body's external linear damping matrix (i.e. non-radiation damping)
    dampingQuad: array
        6x6 array - body's external quadratic damping matrix (i.e. non-radiation damping)
    kHydro: array
        6x6 array - body's hydrostatic stiffness matrix
    kExt: array
        6x6 array - body's additional stiffness matrix

    Returns
    -------
    None

    Raises
    ------
    ValueError
        if no projectDir is passed

    Notes
    -----
    The Hydrostatic.in file is required by HAMS to run, but it may just contain
    zeros if only the hydrodynamic coefficients are of interest.
    '''
    if projectDir is None:
        raise ValueError('No directory has been passed for where to write Hydrostatic.in')
    else:
        projectDir = osp.normpath(osp.join(projectDir, 'Input'))

    f = open(osp.join(projectDir, 'Hydrostatic.in'), 'w')
    f.write(' Center of Gravity:\n ')
    f.write(f'  {cog[0]:10.15E}  {cog[0]:10.15E}  {cog[0]:10.15E} \n')
    f.write(' Body Mass Matrix:\n')
    for i in range(6):
        for j in range(6):
            f.write((f'   {mass[i,j]:10.5E}'))
        f.write('\n')
    f.write(' External Linear Damping Matrix:\n')
    for i in range(6):
        for j in range(6):
            f.write((f'   {dampingLin[i,j]:10.5E}'))
        f.write('\n')
    f.write(' External Quadratic Damping Matrix:\n')
    for i in range(6):
        for j in range(6):
            f.write((f'   {dampingQuad[i,j]:10.5E}'))
        f.write('\n')
    f.write(' Hydrostatic Restoring Matrix:\n')
    for i in range(6):
        for j in range(6):
            f.write((f'   {kHydro[i,j]:10.5E}'))
        f.write('\n')
    f.write(' External Restoring Matrix:\n')
    for i in range(6):
        for j in range(6):
            f.write((f'   {kExt[i,j]:10.5E}'))
        f.write('\n')
    f.close()

def write_control_file(projectDir=None, waterDepth=50.0, incFLim=1, iFType=3, oFType=3, numFreqs=-300,
                       minFreq=0.02, dFreq=0.02, freqList=None, numHeadings=1,
                       minHeading=0.0, dHeading=0.0,
                       refBodyCenter=[0.0, 0.0, 0.0], refBodyLen=1.0, irr=0,
                       numThreads=8):
    '''
    writes HAMS ControlFile.in file

    Parameters
    ----------
    projectDir: str
        main HAMS project directory - ControlFile.in will be written in ./Input/
    waterDepth : float
        water depth (m)
    incFLim: int
        include 0 and infinite frequency limits
            0: do not include
            1: include
    iFType: int
        input frequency type
            1: deepwater wavenumber
            2: finite-depth wavenumber
            3: wave frequency
            4: wave period
            5: wavelength
    oFType: int
        output frequency type
            1: deepwater wavenumber
            2: finite-depth wavenumber
            3: wave frequency
            4: wave period
            5: wavelength
    numFreqs: int
        number of frequencies
            -ve: define minimum frequency and frequency step
            +ve: provide a list of frequencies
    minFreq: float
        minimum frequency (required if numFreqs is -ve)
    dFreq: float
        frequency step (required if numFreqs is -ve)
    freqList: list
        list of frequencies
    numHeadings: int
        number of headings
    minHeading: float
        minimum heading value (degs)
    dHeading: float
        heading step (degs)
    refBodyCenter: list
        reference body center
    refBodyLen: float
        reference body length
    irr: int
        irregular frequency removal option
    numThreads: int
        number of threads

    Returns
    -------

    Raises
    ------

    Notes
    -----

    '''
    if projectDir is None:
        raise ValueError('No directory has been passed for where to write Hydrostatic.in')
    else:
        projectDir = osp.normpath(osp.join(projectDir, './Input'))

    oFileName = 'ControlFile.in'

    f = open(osp.join(projectDir, oFileName), 'w')
    f.write('   --------------HAMS Control file---------------\n\n')
    f.write(f'   Waterdepth  {-waterDepth}D0\n\n')             # note: HAMS expects the depth to be given as negative
    f.write('   #Start Definition of Wave Frequencies\n')
    f.write(f'    0_inf_frequency_limits  {incFLim}\n')
    f.write(f'    Input_frequency_type    {iFType}\n')
    f.write(f'    Output_frequency_type   {oFType}\n')
    f.write(f'    Number_of_frequencies   {numFreqs}\n') # -ve for min & step, +ve for list of frequencies (or periods)
    if not freqList is None and numFreqs > 0:
        f.write(str(np.array(freqList)).replace('[','').replace(']','').replace('\n','')+'\n') # -ve for min & step, +ve for list of frequencies (or periods)
    else:
        f.write(f'    Minimum_frequency_Wmin  {minFreq}D0\n')
        f.write(f'    Frequency_step          {dFreq}D0\n')
    f.write('   #End Definition of Wave Frequencies\n\n')
    f.write('   #Start Definition of Wave Headings\n')
    f.write(f'    Number_of_headings      -{numHeadings}\n')
    f.write(f'    Minimum_heading         {minHeading}D0\n')
    f.write(f'    Heading_step            {dHeading}D0\n')
    f.write('   #End Definition of Wave Headings\n\n')
    f.write(f'    Reference_body_center   {refBodyCenter[0]:.3f} {refBodyCenter[1]:.3f} {refBodyCenter[2]:.3f}\n')
    f.write(f'    Reference_body_length   {refBodyLen}D0\n')
    f.write('    Wave-diffrac-solution   2\n')
    f.write(f'    If_remove_irr_freq      {irr}\n')
    f.write(f'    Number of threads       {numThreads}\n\n')
    f.write('   #Start Definition of Pressure and/or Elevation\n')
    f.write('    Number_of_field_points  0 \n')
    f.write('   #End Definition of Pressure and/or Elevation\n\n')
    f.write('   ----------End HAMS Control file---------------\n')
    f.close()


def read_wamit1(pathWamit1):
    '''
    Read added mass and damping from .1 file (WAMIT format)

    Parameters
    ----------
    pathWamit1 : f-str
        path to .1 file

    Returns
    -------
    addedMass : array
        added mass coefficients (nw x ndof x ndof)
    damping : array
        damping coefficients (nw x ndof x ndof)

    Raises
    ------

    Notes
    -----
    '''
    pathWamit1 = osp.normpath(pathWamit1)
    wamit1 = np.loadtxt(pathWamit1)
    w = np.unique(wamit1[:,0])
    addedMassCol = wamit1[:,3]
    dampingCol = wamit1[:,4]
    addedMass = addedMassCol.reshape((len(w)), 6, 6).transpose(1,2,0)
    damping = dampingCol.reshape((len(w), 6, 6)).transpose(1,2,0)

    return addedMass, damping

def read_wamit1B(pathWamit1, TFlag=0):
    '''
    Read added mass and damping from .1 file (WAMIT format)
    '''
    pathWamit1 = osp.normpath(pathWamit1)
    #wamit1 = np.loadtxt(pathWamit1)
    
    # If they are in the .1 wamit file,
    # Check to make sure that the zero- and infinite-frequency added mass terms are accounted for
    f = open(pathWamit1, 'r')
    # Save the A0 and Ainf data and delete them from the .1 file so that the rest can be read easily
    # Initialize some lists
    freq = []; i = []; j = []; ainf = []; freq0 = []; i0 = []; j0 = []; a0 = [];
    for line in f:                                  # for each line in the .1 file
        if line=="\n":
            pass                                    # skip if the line has nothing in it
        else:
            data = line.split()
            if len(data) < 5:                       # if the line has less than 5 terms, it's either Ainf or A0
                if float(data[0]) < 0.0:            # if the line is an Ainf term
                    freq.append(float(data[0]))
                    i.append(int(data[1]))
                    j.append(int(data[2]))
                    ainf.append(float(data[3]))
                elif float(data[0]) == 0.0:         # if the line is an A0 term
                    freq0.append(float(data[0]))
                    i0.append(int(data[1]))
                    j0.append(int(data[2]))
                    a0.append(float(data[3]))
    # Create the zero- and infinite-frequency added mass matrices
    Ainf = np.zeros([6,6])
    A0 = np.zeros([6,6])
    if ainf:
        for k in range(len(i)):
            Ainf[i[k]-1,j[k]-1] = ainf[k]
        for k in range(len(i0)):
            A0[i0[k]-1,j0[k]-1] = a0[k]
    f.close()                                       # need to close out of the file
    
    f = open(pathWamit1, 'r')                       # and then reopen to allow readlines() to work
    lines = f.readlines()                           # save the whole file into the 'lines' variable
    f.close()                                       # close the file
    
    # count the number of lines that include an A0 or Ainf term
    count = 0
    for line in lines:
        if len(line.split()) < 5:
            count += 1
    
    # delete the lines that include an A0 or Ainf term
    for i in range(count):
        del lines[0]
    

    # reload the rest of the lines variable into a readable matrix for the rest of the function    
    #wamit1 = np.loadtxt(pathWamit1)
    wamit1 = np.loadtxt(lines)
    
    # ------------------------------------------------------------
    
    T = np.unique(wamit1[:,0])
    if TFlag:   # if TFlag=1, the first column values are periods
        w = 2*np.pi/T
    else:       # if TFlag=0, the first column values are frequencies
        w = T
    addedMassCol = wamit1[:,3]
    dampingCol = wamit1[:,4]
    matRow = wamit1[:,1]
    matCol = wamit1[:,2]
    addedMass = np.zeros([len(w),6,6])
    damping = np.zeros([len(w),6,6])
    
    nw = len(w)
    nA = len(addedMassCol)
    x = int(nA/nw)
    for j in range(len(w)):
        for i in range(x):
            addedMass[j,int(matRow[i])-1,int(matCol[i])-1] = addedMassCol[i+x*j]
            damping[j,int(matRow[i])-1,int(matCol[i])-1] = dampingCol[i+x*j]
    
    # transpose(a,b,c) function switches the sizes in an order of k[a] to the 0 spot, k[b] to the 1 spot, and k[c] to the 2 spot
    A = addedMass.transpose(1,2,0) # transposes addedMass.shape = (100,6row,6col) to (6row,6col,100)
    B = damping.transpose(1,2,0)
    
    return A, B, w      # can return Ainf and A0 here too

def read_wamit3(pathWamit3):
    '''
    Read excitation force coefficients from .3 file (WAMIT format)

    Parameters
    ----------
    pathWamit3 : f-str
        path to .3 file

    Returns
    -------
    mod:

    phase:

    real:

    imag:

    Raises
    ------

    Notes
    -----
    '''
    pathWamit3 = osp.normpath(pathWamit3)
    wamit3 = np.loadtxt(pathWamit3)
    w = np.unique(wamit3[:,0])
    headings = np.unique(wamit3[:,1])
    mod = wamit3[:,3].reshape((len(w), 6)).transpose(1,0)
    phase = wamit3[:,4].reshape((len(w), 6)).transpose(1,0)
    real = wamit3[:,5].reshape((len(w), 6)).transpose(1,0)
    imag = wamit3[:,6].reshape((len(w), 6)).transpose(1,0)

    return mod, phase, real, imag

def read_wamit3B(pathWamit3, TFlag=0):
    '''
    Read excitation force coefficients from .3 file (WAMIT format)
    '''
    
    pathWamit3 = osp.normpath(pathWamit3)
    wamit3 = np.loadtxt(pathWamit3)
    
    if TFlag:   # if TFlag=1, the first column values are periods
        T = np.unique(wamit3[:,0])
        w = 2*np.pi/np.flip(T)
    else:       # if TFlag=0, the first column values are frequencies
        w = np.unique(wamit3[:,0])
        
    headings = np.unique(wamit3[:,1])
        
    modCol = wamit3[:,3]
    phaseCol = wamit3[:,4]
    realCol = wamit3[:,5]
    imagCol = wamit3[:,6]
    
    mod = np.zeros([len(w), len(headings), 6])
    phase = np.zeros([len(w), len(headings), 6])
    real = np.zeros([len(w), len(headings), 6])
    imag = np.zeros([len(w), len(headings), 6])
    
    nw = len(w)
    nh = len(headings)
    nF = len(modCol)
    x = int(nF/nh/nw)
    y = int(nF/nw)
    
    for i in range(len(w)):
        for j in range(len(headings)):
            for k in range(x):
                mod[i,j,k] = modCol[k+x*j+y*i]
                phase[i,j,k] = phaseCol[k+x*j+y*i]
                real[i,j,k] = realCol[k+x*j+y*i]
                imag[i,j,k] = imagCol[k+x*j+y*i]
    
    M = mod.transpose(1,2,0) # transposes a mod matrix of size (100, 37, 6) to (37, 6, 100)
    P = phase.transpose(1,2,0)
    R = real.transpose(1,2,0)
    I = imag.transpose(1,2,0)
    
    return M, P, R, I, w, headings

def run_hams(projectDir):
    '''call the HAMS_x64.exe program in the specified project directory'''
    # get absolute path to the local HAMS_x64.exe program
    hamsDir = osp.dirname(__file__)
    hamsExe = './bin/HAMS_x64.exe'
    hamsPath = osp.abspath(osp.join(hamsDir, hamsExe))
    # change directory to where the HAMS input files are
    workingDir = os.getcwd()
    os.chdir(projectDir)
    # run HAMS
    if platform.system() == 'Windows':
        sub.run([f'{hamsPath}'])
    else:
        hams.hams_lib.exec()
    # change back to working directory
    os.chdir(workingDir)
