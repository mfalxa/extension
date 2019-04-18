import yaml
import io


def genConfigFileLTCube(name, eventFilePath, spacecraftFilePath, ltcubeFilePath, targetGLon, targetGLat, roiwidth, binsz, binsperdec, emin, emax, src_roiwidth) :

	config = {
		'data' : {
			'evfile' : eventFilePath,
			'scfile' : spacecraftFilePath,
			'ltcube' : ltcubeFilePath
		},
		'binning' : {
			'coordsys' : 'GAL',
			'roiwidth' : roiwidth,
			'binsz' : binsz,
			'binsperdec' : binsperdec
		},
		'selection' : {
			'emin' : emin,
			'emax' : emax,
			'zmax' : 90,
			'evclass' : 128,
			'evtype' : 3,
			'filter' : None,
			'glon' : targetGLon,
			'glat' : targetGLat
		},
		'gtlike' : {
			'edisp' : True,
			'irfs' : 'P8R3_SOURCE_V2'
		},
		'model' : {
			'src_roiwidth' : src_roiwidth,
			'galdiff' : '$FERMI_DIFFUSE_DIR/gll_iem_v07.fits',
			'isodiff': 'iso_P8R3_SOURCE_V2_v1.txt',
			'catalogs' : ['4FGL']
		},
		'fileio' : {
			'outdir' : '/users-data/mfalxa/code/' + name
		}
	}

	with io.open('/users-data/mfalxa/code/' + name + '/config.yaml', 'w') as outfile:
		yaml.dump(config, outfile, default_flow_style=False)


def genConfigFile(name, eventFilePath, spacecraftFilePath, targetGLon, targetGLat, roiwidth, binsz, binsperdec, emin, emax, src_roiwidth) :

	config = {
		'data' : {
			'evfile' : eventFilePath,
			'scfile' : spacecraftFilePath
		},
		'binning' : {
			'coordsys' : 'GAL',
			'roiwidth' : roiwidth,
			'binsz' : binsz,
			'binsperdec' : binsperdec
		},
		'selection' : {
			'emin' : emin,
			'emax' : emax,
			'zmax' : 90,
			'evclass' : 128,
			'evtype' : 3,
			'filter' : None,
			'glon' : targetGLon,
			'glat' : targetGLat
		},
		'gtlike' : {
			'edisp' : True,
			'irfs' : 'P8R3_SOURCE_V2'
		},
		'model' : {
			'src_roiwidth' : src_roiwidth,
			'galdiff' : '$FERMI_DIFFUSE_DIR/gll_iem_v07.fits',
			'isodiff': 'iso_P8R3_SOURCE_V2_v1.txt',
			'catalogs' : ['4FGL']
		},
		'fileio' : {
			'outdir' : '/users-data/mfalxa/code/' + name
		}
	}

	with io.open('/users-data/mfalxa/code/' + name + '/config.yaml', 'w') as outfile:
		yaml.dump(config, outfile, default_flow_style=False)
