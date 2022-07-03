import os, sys
import numpy as np
from astropy import table
from astropy.io import fits
from astLib import *
from nemo import completeness, plotSettings, catalogs, signals
import pylab as plt
from scipy import stats
import yaml
import glob


def make_truth_mock(mode, configdict):

	if mode == 'with_Q':
		# Make a 'true' mock - use the truth catalog, get true_SNR by looking up noise in the selFn dir
		truthTab=table.Table().read(configdict['path2truthcat'])
		noiseMapFileName= configdict['path2noisemap']
		with fits.open(noiseMapFileName) as img:
		    for ext in img:
		        if ext.data is not None:
		            break
		    rmsMap=ext.data
		    wcs=astWCS.WCS(ext.header, mode = 'pyfits')
		truthTab['true_SNR']=0.0
		truthTab['missed']=0 # For spotting (a handful) of clusters that fell outside mask
		for row in truthTab:
		    if wcs.coordsAreInImage(row['RADeg'], row['decDeg']) is True:
		        x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
		        x=int(round(x)); y=int(round(y))
		        if x < rmsMap.shape[1]-1 and y < rmsMap.shape[0]-1 and rmsMap[y, x] != 0:
		            row['true_SNR']=(row['true_fixed_y_c']*1e-4) / rmsMap[y, x]
		        else:
		            row['missed']=1
		    else:
		        row['missed']=1

	elif mode == 'without_Q':	 

		selFn=completeness.SelFn(configdict['path2selFn'], SNRCut = configdict['predSNRCut'], zStep = configdict['selFnZStep'],
                     enableDrawSample = configdict['makeMock'], massFunction = configdict['massFunc'],
                     applyRelativisticCorrection = configdict['relativisticCorrection'],
                     rhoType = configdict['rhoType'], delta = configdict['delta'])       

		truthTab=table.Table().read(configdict['path2truthcat'])
		noiseMapFileName = configdict['path2noisemap']

		with fits.open(noiseMapFileName) as img:
		    for ext in img:
		        if ext.data is not None:
		            break
		    rmsMap=ext.data
		    wcs=astWCS.WCS(ext.header, mode = 'pyfits')
		Q=signals.QFit(configdict['path2Qfunc'])
		truthTab['true_SNR']=0.0
		truthTab['true_fixed_y_c']=0.0
		truthTab['true_Q']=0.0
		truthTab['missed']=0 # For spotting (a handful) of clusters that fell outside mask
		print("WARNING: We don't have true_fixed_y_c or true_Q - we reconstruct those here.")
		for row in truthTab:
		    if wcs.coordsAreInImage(row['RADeg'], row['decDeg']) is True:
		        x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
		        x=int(round(x)); y=int(round(y))
		        if x < rmsMap.shape[1]-1 and y < rmsMap.shape[0]-1 and rmsMap[y, x] != 0:
		            # Need to know tileNames for objects in truth catalog
		            thisQ=Q.getQ(signals.calcTheta500Arcmin(row['redshift'], row['true_M500c']*1e14,
		                                                    selFn.mockSurvey.cosmoModel), tileName = row['tileName'])
		            #Ez=ccl.h_over_h0(signals.fiducialCosmoModel, 1/(1+row['redshift']))
		            row['true_Q']=thisQ
		            row['true_fixed_y_c']=row['true_y_c']*thisQ
		            row['true_SNR']=(row['true_fixed_y_c']*1e-4) / rmsMap[y, x]
		        else:
		            row['missed']=1
		    else:
		        row['missed']=1

	return truthTab


def make_nemo_mock(configdict):

	selFn=completeness.SelFn(configdict['path2selFn'], SNRCut = configdict['predSNRCut'], zStep = configdict['selFnZStep'],
                         enableDrawSample = configdict['makeMock'], massFunction = configdict['massFunc'],
                         applyRelativisticCorrection = configdict['relativisticCorrection'],
                         rhoType = configdict['rhoType'], delta = configdict['delta'])

	mockTab=selFn.generateMockSample(mockOversampleFactor = configdict['predAreaScale'],
                                     applyPoissonScatter = configdict['applyPoissonScatter'])

	return mockTab


def get_nemo_pred(configdict, zbins):

	selFn=completeness.SelFn(configdict['path2selFn'], SNRCut = configdict['predSNRCut'], zStep = configdict['selFnZStep'],
                     enableDrawSample = configdict['makeMock'], massFunction = configdict['massFunc'],
                     applyRelativisticCorrection = configdict['relativisticCorrection'],
                     rhoType = configdict['rhoType'], delta = configdict['delta'])

	predMz=selFn.compMz*selFn.mockSurvey.clusterCount
	#TODO: Ask Matt where the minimal mass comes from
	predNz_fineBins=predMz[:, np.greater(selFn.mockSurvey.log10M, np.log10(5e13))].sum(axis = 1)

	predNz=np.zeros(zbins.shape[0]-1)
	for i in range(len(zbins)-1):
	    zMin=zbins[i]
	    zMax=zbins[i+1]
	    mask=np.logical_and(selFn.mockSurvey.z > zMin, selFn.mockSurvey.z <= zMax)
	    predNz[i]=predNz_fineBins[mask].sum()

	return predNz


def bin_catalog(catalog, zbins, qbins, SNR_tag='SNR'):

	# redshift bins for N(z)
	zarr = 0.5*(zbins[:-1] + zbins[1:])
	qarr = 0.5*(qbins[:1] + qbins[1:])

	delN2Dcat, _, _ = np.histogram2d(catalog['redshift'], catalog[SNR_tag], bins=[zbins, qbins])

	return delN2Dcat, zarr, qarr




