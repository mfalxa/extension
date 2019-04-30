import matplotlib

matplotlib.use('Agg')

import numpy as np
import yaml
from fermipy.gtanalysis import GTAnalysis
import genConfig
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
 
def distanceFromCenter(centerGLon, centerGLat, sourceObject) :
	center = SkyCoord(centerGLon * u.degree, centerGLat * u.degree)
	source = SkyCoord(sourceObject['glon'] * u.degree, sourceObject['glat'] * u.degree)
	separation  = center.separation(source)
	return separation.value

class ExtensionFit :

	def __init__(self, configFile) :
		
		self.gta = GTAnalysis(configFile, logging={'verbosity' : 3})
                self.target = None

	''' INITIALIZE '''

	def initialize(self, sizeROI, debug, addToROI) :

		self.gta.setup()
                if self.gta.config['selection']['emin'] >= 10000 :
                        self.gta.set_parameter('galdiff', 'Scale', 30000)
                                      
		# Get model source names
		sourceList = self.gta.get_sources(exclude=['isodiff', 'galdiff'])

		catalog = fits.get_data('/users-data/mfalxa/code/gll_psch_v13.fit', 1)

		# Delete sources unassociated with TS < 50
		for i in range(len(sourceList)) :
			if sourceList[i]['catalog']['TS_value'] < 50. and catalog['CLASS'][catalog['Source_Name'] == sourceList[i]['name']][0] == '' :
				self.gta.delete_source(sourceList[i]['name'])

                closests = self.gta.get_sources(distance=rInner, exclude=['isodiff', 'galdiff'])

		# Delete all unidentified sources
		for i in range(len(closests)) :
			if catalog['CLASS'][catalog['Source_Name'] == closests[i]['name']][0].isupper() == False  :
				self.gta.delete_source(closests[i]['name'])
			if catalog['CLASS'][catalog['Source_Name'] == closests[i]['name']][0] == 'SFR' :
				self.target = closests[i]['name']

                # If debug, save ROI and make plots
		if debug == True :
			self.gta.write_roi('startModel')
			self.gta.residmap(prefix='start', make_plots=True)
			self.gta.make_plots('start')

		# Optmize spectral parameters for sources with npred > 1
		self.gta.optimize(npred_threshold=1, skip=['isodiff'])

		# Get model source names
		sourceList = self.gta.get_sources(distance=sizeROI + addToROI, square=True, exclude=['isodiff', 'galdiff'])

		# Iterate source localizing on source list
		for i in range(len(sourceList)) :
			if sourceList[i].extended == False :
				self.gta.localize(sourceList[i]['name'], write_fits=False, write_npy=False, update=True)

		# Free sources within ROI size + extra distance from center
		self.gta.free_sources(distance=sizeROI + addToROI, square=True)

		# Re-optimize ROI
		self.gta.optimize(skip=['isodiff'])

		# Save and make plots if debug
		if debug == True :
			self.gta.write_roi('modelInitialized')
			self.gta.residmap(prefix='initialized', make_plots=True)
			self.gta.make_plots('initialized')

		# Lock sources
		self.gta.free_sources(free=False)

	''' OUTER REGION '''

	def outerRegionAnalysis(self, centerGLon, centerGLat, sizeROI, rInner, sqrtTsThreshold, minSeparation, debug) :

		self.gta.free_sources(distance=sizeROI, pars='norm', square=True, free=True)
		self.gta.free_sources(distance=rInner, free=False)
		self.gta.free_source('galdiff', free=True)
		self.gta.free_source('isodiff', free=False)

		# Seek new sources until none are found
		sourceModel = {'SpectrumType' : 'PowerLaw', 'Index' : 2.0, 'Scale' : 30000, 'Prefactor' : 1.e-15, 'SpatialModel' : 'PointSource'}
		newSources = self.gta.find_sources(sqrt_ts_threshold=sqrtTsThreshold, min_separation=minSeparation, model=sourceModel,
				**{'search_skydir' : self.gta.roi.skydir, 'search_minmax_radius' : [rInner, sizeROI]})
			
		if len(newSources) > 0 :
			for i in range(len(newSources)) :
				if newSources['sources'][i]['ts'] > 100. :
					self.gta.set_source_spectrum(newSources['sources'][i]['name'], spectrum_type='LogParabola')
					self.gta.free_source(newSources['sources'][i]['name'])
					self.gta.fit()
                                        self.gta.free_sources(newSources['sources'][i]['name'], free=False) 

		# Optimize all ROI
		self.gta.optimize(skip=['isodiff'])

		# Save sources found
		if debug == True :
			np.save('sourcesFoundOuter.npy', sourcesFoundDict)
			self.gta.residmap(prefix='outer', make_plots=True)
			self.gta.write_roi('outerAnalysisROI')
			self.gta.make_plots('outer')


	''' INNER REGION '''

	def innerRegionAnalysis(self, sizeROI, rInner, maxIter, sqrtTsThreshold, minSeparation, debug) :

		self.gta.free_sources(distance=sizeROI, square=True, free=False)
		self.gta.free_sources(distance=rInner, free=True, exclude=['isodiff'])

		# Keep closest source if identified with star forming region in catalog or look for new source closest to center within Rinner
		if self.target != None :
			print('Closest source identified with star forming region : ', self.target)
			self.gta.set_source_morphology(self.target, **{'spatial_model' : 'PointSource'})
		else :
			closeSources = self.gta.find_sources(sqrt_ts_threshold=2., min_separation=minSeparation, max_iter=1,
							**{'search_skydir' : self.gta.roi.skydir, 'search_minmax_radius' : [0., rInner]})
			dCenter = np.array([])
			for i in range(len(closeSources['sources'])) :
				dCenter = np.append(dCenter, self.gta.roi.skydir.separation(closeSources['sources'][i].skydir).value)
			self.target = closeSources['sources'][np.argmin(dCenter)]['name']
			for i in [x for x in range(len(closeSources['sources'])) if x != (np.argmin(dCenter))] :
				self.gta.delete_source(closeSources['sources'][i]['name'])
			self.gta.optimize(skip=['isodiff'])

		# Initialize n sources array
		nSources = []

                # Save ROI without extension fit
		self.gta.write_roi('nSourcesFit')
		
		if debug == True :
			self.gta.make_plots('innerInit')
                        self.gta.residmap(prefix='innerInit', make_plots=True)
		
		# Test for extension
		extensionTest = self.gta.extension(self.target, make_plots=True, write_npy=debug, write_fits=debug, spatial_model='RadialDisk', update=True, free_background=True, fit_position=True)
		extLike = extensionTest['loglike_ext']
		extAIC = 2 * (len(self.gta.get_free_param_vector()) - extLike)
		self.gta.write_roi('extFit')

		if debug == True :
			self.gta.make_plots('ext0')

		self.gta.load_roi('nSourcesFit', reload_sources=True)	
	
		for i in range(1, maxIter + 1) :
			
			# Test for n point sources
			nSourcesTest = self.gta.find_sources(sources_per_iter=1, sqrt_ts_threshold=sqrtTsThreshold, min_separation=minSeparation, max_iter=1,
							**{'search_skydir' : self.gta.roi.skydir, 'search_minmax_radius' : [0., rInner]})

			if len(nSourcesTest['sources']) > 0 :

				if nSourcesTest['sources'][0]['ts'] > 100. :
					self.gta.set_source_spectrum(nSourcesTest['sources'][0]['name'], spectrum_type='LogParabola')
					self.gta.free_source(nSourcesTest['sources'][0]['name'])
					self.gta.fit()
                                        self.gta.free_sources(nSourcesTest['sources'][0]['name'], free=False)

				if debug == True :
					self.gta.make_plots('nSources' + str(i))

				nSources.append(nSourcesTest['sources'])
				self.gta.localize(nSourcesTest['sources'][0]['name'], write_npy=debug, write_fits=debug, update=True)
				nAIC = 2 * (len(self.gta.get_free_param_vector()) - self.gta._roi_data['loglike'])
                                self.gta.free_source(nSourcesTest['sources'][0]['name'])
				self.gta.write_roi('nSourcesFit')
				
				# Estimate Akaike Information Criterion difference between both models
				dm = extAIC - nAIC
				print('AIC difference between both models = ', dm)

				# Estimate TS_m+1
				extensionTestPlus = self.gta.extension(self.target, make_plots=True, write_npy=debug, write_fits=debug, spatial_model='RadialDisk', update=True, free_background=True, fit_position=True)
				TSm1 = 2 * (extensionTestPlus['loglike_ext'] - extLike)
				print('TSm+1 = ', TSm1)

				if debug == True :
					self.gta.make_plots('ext' + str(i))

				if dm < 0 and TSm1 < 16 :
					self.gta.load_roi('extFit', reload_sources=True)
					break
				else :

					# Set extension test to current state and save current extension fit ROI and load previous nSources fit ROI
					extentionTest = extensionTestPlus
					extLike = extensionTestPlus['loglike_ext']
					extAIC = 2 * (len(self.gta.get_free_param_vector()) - extLike)
					self.gta.write_roi('extFit')
					self.gta.load_roi('nSourcesFit', reload_sources=True)
			
			else :
				self.gta.load_roi('extFit', reload_sources=True)
				break

		self.gta.fit()
                
	''' CHECK OVERLAP '''
