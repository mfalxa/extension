import matplotlib

matplotlib.use('Agg')

import numpy as np
import os
import yaml
import genConfig
import statistictica
import astropy.io.fits as fits
from fermipy.gtanalysis import GTAnalysis


''' GLOBAL VARIABLES '''

debug = True
newListFile = False

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
	for j in range(len(clusters[i]) - 3) :


		''' SET PARAMETERS '''
		
		# Config file
		name = clusters[i][j]['name']
		radius = clusters[i][j]['radius']
		targetGLon = clusters[i][j]['galLon']
		targetGLat = clusters[i][j]['galLat']
		roiwidth = 6. + 2 * radius
		binsz = 0.1
		binsperdec = 25
		emin = 10000
		emax = 1000000
		src_roiwidth = 10. + 2 * radius

		# Analysis
		sizeROI = roiwidth / 2
		addToROI = 0.1
		rInner = radius + 1.
		sqrtTsThreshold = 3.
		maxIter = 5
		minSeparation = 0.5

		# TS conditions		
		TSMin = 50.
		dmMin = 0.
		TSm1Min = 16.

		# SED
		logeBins = np.linspace(np.log10(emin), np.log10(emax), 11)

		if radius > 1. :
			pass


		''' WRITE VARIABLES IN *.TXT FILE '''

		parameters = np.array(['name : '+str(name), 'glon : '+str(targetGLon), 'glat : '+str(targetGLat), 'roiwidth : '+str(roiwidth), 'binsz : '+str(binsz), 'binsperdec : '+str(binsperdec), 'emin : '+str(emin), 'emax : '+str(emax), 'src_roiwidth : '+str(src_roiwidth), 'sizeROI : '+str(sizeROI), 'addToROI : '+str(addToROI), 'radius : '+str(radius), 'rInner : '+str(rInner), 'sqrtTsThreshold : '+str(sqrtTsThreshold), 'maxIter : '+str(maxIter), 'minSeparation : '+str(minSeparation), 'TSmin : '+str(TSMin), 'dmMin : '+str(dmMin), 'TSm1Min : '+str(TSm1Min)])
		print(parameters)
		np.savetxt('/users-data/mfalxa/code/' + name + '/parameters.txt', parameters, fmt='%s')



		''' GENERATE CONFIG FILE AND GTA INSTANCE '''
		
		os.system('mkdir ' + name)
		
		if newListFile == False :

			genConfig.genConfigFileLTCube(name,
							'/users-data/mfalxa/code/Cyg_2704_07h02/ft1_00.fits',
							'/users-data/mfalxa/code/Cyg_2704_07h02/scfile_00.txt',
							'/users-data/mfalxa/code/Cyg_2704_07h02/ltcube_00.fits',
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
	
		Analysis = statistictica.ExtensionFit('/users-data/mfalxa/code/' + name + '/config.yaml')


#==================================================#
		''' RUN ANALYSIS '''
#==================================================#

		Analysis.initialize(sizeROI, rInner, addToROI, TSMin, debug)
		Analysis.outerRegionAnalysis(sizeROI, rInner, sqrtTsThreshold, minSeparation, debug)
		Analysis.innerRegionAnalysis(sizeROI, rInner, maxIter, sqrtTsThreshold, minSeparation, dmMin, TSm1Min, debug)
		Analysis.overlapDisk(radius)
		Analysis.gta.write_roi(name)
		Analysis.gta.make_plots('end')
		Analysis.gta.sed(name=Analysis.target['name'], loge_bins=logeBins, make_plots=True)
		Analysis.gta.residmap(prefix='end', make_plots=True)
