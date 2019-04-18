import numpy as np
import yaml
from fermipy.gtanalysis import GTAnalysis
import genConfig
from astropy import units as u
from astropy.coordinates import SkyCoord

 
def distanceFromCenter(centerGLon, centerGLat, sourceObject) :
	center = SkyCoord(centerGLon * u.degree, centerGLat * u.degree)
	source = SkyCoord(sourceObject['glon'] * u.degree, sourceObject['glat'] * u.degree)
	separation  = center.separation(source)
	return separation.value

class ExtensionFit :

	def __init__(self, configFile) :
		
		self.gta = GTAnalysis(configFile, logging={'verbosity' : 3})		

	''' INITIALIZE '''

	def initialize(self, sizeROI, debug, addToROI) :

		self.gta.setup()

		# If debug, save ROI and make plots
		if debug == True :
			self.gta.write_roi('startModel')
			self.gta.make_plots('start')

		# Get model source names
		sourceList = self.gta.get_sources(distance=sizeROI + addToROI, square=True, exclude=['isodiff', 'galdiff'])

		# Delete sources unassociated with TS < 100 or flags 5, 6, 8
		for i in range(len(sourceList)) :
			if sourceList[i]['catalog']['Flags'] in (5, 6, 8) :
				self.gta.delete_source(sourceList[i]['name'])
			if sourceList[i]['catalog']['CLASS'].isupper() == False and sourceList[i]['catalog']['TS'] < 100 :
				self.gta.delete_source(sourceList[i]['name'])

		# Optmize spectral parameters for sources with npred > 1
		self.gta.optimize(npred_threshold=1)

		# Get model source names
		sourceList = self.gta.get_sources(distance=sizeROI + addToROI, square=True, exclude=['isodiff', 'galdiff'])

		# Iterate source localizing on source list
		for i in range(len(sourceList)) :
			if sourceList[i].extended == False :
				self.gta.localize(sourceList[i]['name'], write_fits=False, write_npy=False, update=True)

		# Free sources within ROI size + extra distance from center
		self.gta.free_sources(distance=sizeROI + addToROI, square=True)

		# Re-optimize ROI
		self.gta.optimize()

		# Lock sources
		self.gta.free_sources(distance=sizeROI + addToROI, square=True, free=False)

		if debug == True :
			self.gta.write_roi('modelInitialized')
			self.gta.make_plots('initialized')

	''' OUTER REGION '''

	def outerRegionAnalysis(self, centerGLon, centerGLat, sizeROI, rInner, sqrtTsThreshold, minSeparation, debug) :

		self.gta.free_sources(distance=sizeROI, square=True, free=False)

		# Seek new sources until none are found
		sourcesFound = 1
		sourcesFoundDict = np.array([])
		while sourcesFound != 0 :
			newSources = self.gta.find_sources(sqrt_ts_threshold=sqrtTsThreshold, min_separation=minSeparation, max_iter=1,
				**{'search_skydir' : self.gta.roi.skydir, 'search_minmax_radius' : [rInner, sizeROI]})
			sourcesFoundDict = np.append(sourcesFoundDict, newSources['sources'])
			sourcesFound = len(newSources['sources'])

		# Optimize all ROI
		self.gta.optimize()

		# Save sources found
		if debug == True :
			np.save('sourcesFoundOuter.npy', sourcesFoundDict)
			self.gta.write_roi('outerAnalysisROI')
			self.gta.make_plots('outer')


	''' INNER REGION '''

	def innerRegionAnalysis(self, sizeROI, rInner, maxIter, sqrtTsThreshold, minSeparation, debug) :

		self.gta.free_sources(distance=sizeROI, square=True, free=False)
		
		# Find closest sources to ROI center
		closest = None
		closests = self.gta.get_sources(distance=rInner, exclude=['isodiff', 'galdiff'])

		# Delete all unidentified sources
		for i in range(len(closests)) :
			if closests[i]['catalog']['CLASS'].isupper() == False :
				self.gta.delete_source(closests[i]['name'])
			if closests[i]['catalog']['CLASS'] == 'SFR' :
				closest = closests[i]['name']
		
		# Keep closest source if identified with star forming region in catalog or look for new source closest to center within Rinner
		if closest != None :
			print('Closest source identified with star forming region : ', closest)
			self.gta.set_source_morphology(closest, **{'spatial_model' : 'PointSource'})
			self.gta.optimize()
		else :
			closeSources = self.gta.find_sources(sqrt_ts_threshold=2., min_separation=minSeparation, max_iter=1,
							**{'search_skydir' : self.gta.roi.skydir, 'search_minmax_radius' : [0., rInner]})
			dCenter = np.array([])
			center = self.gta.roi.skydir
			for i in range(len(closeSources['sources'])) :
				dCenter = np.append(dCenter, center.separation(closeSources['sources'][i].skydir).value)
			closest = closeSources['sources'][np.argmin(dCenter)]['name']
			for i in [x for x in range(len(closeSources['sources'])) if x != (np.argmin(dCenter))] :
				self.gta.delete_source(closeSources['sources'][i]['name'])
			self.gta.optimize()

		# Initialize n sources array
		nSources = np.array([])

                # Save ROI without extension fit
		self.gta.write_roi('nSourcesFit')
		
		if debug == True :
			self.gta.make_plots('innerInit')
		
		# Test for extension
		extensionTest = self.gta.extension(closest, write_npy=debug, write_fits=debug, free_background=False, spatial_model='RadialDisk', update=True)
		extLike = self.gta.like()
		extAIC = 2 * (len(self.gta.get_free_param_vector()) - extLike)
		self.gta.write_roi('extFit')

		if debug == True :
			self.gta.make_plots('ext0')

		
		for i in range(1, maxIter + 1) :
			
			self.gta.load_roi('nSourcesFit')

			# Test for n point sources
			nSourcesTest = self.gta.find_sources(sources_per_iter=1, sqrt_ts_threshold=sqrtTsThreshold, min_separation=minSeparation, max_iter=1,
							**{'search_skydir' : self.gta.roi.skydir, 'search_minmax_radius' : [0., rInner]})
			
			if len(nSourcesTest['sources']) > 0 : 
			
				if debug == True :
					self.gta.make_plots('nSources' + str(i))

				nSources = np.append(nSources, nSourcesTest['sources'])
				self.gta.free_source(nSourcesTest['sources'][0]['name'])
				nFit = self.gta.fit()
				self.gta.localize(nSourcesTest['sources'][0]['name'], write_npy=debug, write_fits=debug, update=True)
				nAIC = 2 * (len(self.gta.get_free_param_vector()) - self.gta.like())
				self.gta.write_roi('nSourcesFit')
				
				# Estimate Akaike Information Criterion difference between both models
				dm = extAIC - nAIC
				print('AIC difference between both models = ', dm)

				# Estimate TS_m+1
				extensionTestPlus = self.gta.extension(closest, write_npy=debug, write_fits=debug, free_background=False, spatial_model='RadialDisk', update=True)
				TSm1 = 2 * (self.gta.like() - extLike)
				print('TSm+1 = ', TSm1)

				if debug == True :
					self.gta.make_plots('ext' + str(i))

				if dm < 0 and TSm1 < 16 :
					self.gta.load_roi('extFit')
					break
				else :

					# Set extension test to current state and save current extension fit ROI and load previous nSources fit ROI
					extentionTest = extensionTestPlus
					extLike = self.gta.like()
					extAIC = 2 * (len(self.gta.get_free_param_vector()) - extLike)
					self.gta.write_roi('extFit')
					self.gta.load_roi('nSourcesFit')
			
			else :
				self.gta.load_roi('extFit')
				break


	''' CHECK OVERLAP '''
