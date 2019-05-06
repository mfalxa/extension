import numpy as np
import os
import yaml
import genConfig
import statisticalTest
from fermipy.gtanalysis import GTAnalysis


''' GLOBAL VARIABLES '''

debug = True
newListFile = True

''' CREATE *.LST EVENT FILE AND SPACECRAFT FILE '''

if newListFile == True :
	os.system('ls /project-data/glast/flight/fssc/p8v3/ft1/*.fits > ft1.lst')
	os.system('ls /project-data/glast/flight/fssc/p8v3/ft2/*.fits > ft2.lst')

''' LOAD CLUSTER CATALOG '''

clusters = yaml.load(open('clusterz.yaml'))

#########################################
''' ITERATE ANALYSIS FOR EACH CLUSTER '''
#########################################

a = [10]

for i in a :
	for j in range(len(clusters[i])) :


		''' SET PARAMETERS '''
		
		# Config file
		name = clusters[i][j]['name']
		targetGLon = clusters[i][j]['galLon']
		targetGLat = clusters[i][j]['galLat']
		roiwidth = 6.
		binsz = 0.1
		binsperdec = 25
		emin = 10000
		emax = 1000000
		src_roiwidth = 10.

		# Analysis
		sizeROI = roiwidth / 2
		addToROI = 0.1
		radius = clusters[i][j]['radius']
		rInner = radius + 1.
		sqrtTsThreshold = 3.
		maxIter = 5
		minSeparation = 0.5

		if radius > 1. :
			pass

		''' WRITE VARIABLES IN *.TXT FILE '''

		parameters = np.array(['name : '+str(name), 'glon : '+str(targetGLon), 'glat : '+str(targetGLat), 'roiwidth : '+str(roiwidth), 'binsz : '+str(binsz), 'binsperdec : '+str(binsperdec), 'emin : '+str(emin), 'emax : '+str(emax), 'src_roiwidth : '+str(src_roiwidth), 'sizeROI : '+str(sizeROI), 'addToROI : '+str(addToROI), 'radius : '+str(radius), 'rInner : '+str(rInner), 'sqrtTsThreshold : '+str(sqrtTsThreshold), 'maxIter : '+str(maxIter), 'minSeparation : '+str(minSeparation)])
		np.savetxt('parameters.txt', parameters)

		''' GENERATE CONFIG FILE AND GTA INSTANCE '''
		
		os.system('mkdir ' + name)
		
		if newListFile == False :

			genConfig.genConfigFileLTCube(name,
							'ft1.lst',
							'ft2.lst',
							'/users-data/mfalxa/code/' + name + '/ltcube_00.fits',
							targetGLon,
							targetGLat,
							roiwidth,
							binsz,
							binsperdec,
							emin,
							emax,
							src_roiwidth)
		else :
			genConfig.genConfigFile(name,
						'ft1.lst',
						'ft2.lst',
						targetGLon,
						targetGLat,
						roiwidth,
						binsz,
						binsperdec,
						emin,
						emax,
						src_roiwidth)
	
		Analysis = statisticalTest.ExtensionFit('/users-data/mfalxa/code/' + name + '/config.yaml')

#==================================================#
		''' RUN ANALYSIS '''
#==================================================#

		Analysis.initialize(sizeROI, addToROI)
		Analysis.outerRegionAnalysis(targetGLon, targetGLat, sizeROI, rInner, sqrtTsThreshold, minSeparation, debug)
		Analysis.innerRegionAnalysis(sizeROI, rInner, maxIter, sqrtTsThreshold, minSeparation, debug)
		Analysis.gta.write_roi(name)
