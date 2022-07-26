'''
 Calculates the topographic index (a/ tan(beta))
  
 author: Mino Sorribas (mino.sorribas@gmail.com)
 2022.07.20
 
 Pythonized version of original C code from:
 https://github.com/ICHydro/topmodel/blob/master/topmodel/src/c_topidx.c
 
 I recommend the following papers
 
 K. J. BEVEN & M. J. KIRKBY (1979) A physically based, variable
 contributing area model of basin hydrology / Un modèle à base physique de zone d'appel
 variable de l'hydrologie du bassin versant, Hydrological Sciences Journal, 24:1, 43-69
 DOI: 10.1080/02626667909491834
 
 P QUINN P, K BEVEN et al (1991),  The prediction of hillslope flow paths
 for distributed hydrological modelling using digital terrain models.
 Hydrol. Process., 5: 59-79. https://doi.org/10.1002/hyp.3360050106
 
 EM O’LOUGHLIN (1986). Prediction of Surface Saturation Zones in Natural Catchments
 by Topographic Analysis. Water Resources Research, 22(5), 794–804.
 doi:10.1029/wr022i005p00794 

Issues:

    - I do not understand WHY rivers/sink pixels are allowed to be
    processed even if upslopes were not resolved, which is a clear
    restriction for the "upland" pixels.
    
    River/sink pixel has the following code for the topographic index.
        atb[i,j] = np.log(area[i,j] / (2 * sumtb))
    ..thus it seems to me its dependant of proper flow accumulation upstream    
   
 
Notes:
     
    - According to Beven & Kirkby (2009) the 'a' parameter is the
     "area drained per unit contour length", but i'm not sure if people
     elsewhere are doing it properly, so here I am.
     
    - By translating the original code, I understand that the
    (upstream) area is updated during the process resulting in a
    flow_acc which accounts for a weighted contour length adjustment.
    I understand it's is based on Quinn et al. (1991)
            
    - Its not clear to me, but maybe this implementation requires a
    conditioned dem (sinks+fdr), so using raw or bare-earth dem could
    put additional challenge.  i'm too not sure, yet.
  
    - I included lots of comments to make it more understandable (to me)
    
    - Adapted running window indexes and conditions (jj,ii,etc.)

    - Replaced variable ZERO = 0.0000001 to 0.
    - Replaced initial value for atb from -9.9 to -1
    - Final table is filtered for non-negative atb
    - Returns both atb and area

 
    - (not really) issues carried from original code:
      : average downslope slope is also calculated but not returned 
 
 
TODO:
    - assert dem and river matrixes
    - get resolution from raster metadata
    - test! :/
 
'''

import numpy as np
import rasterio
import pandas as pd

# read dem and river arrays
with rasterio.open('dem.tif','r') as src:
    dem = src.read(1)
with rasterio.open('river.tif','r') as src:
    river = src.read(1)


# metadata
nrow = dem.height
ncol = dem.width
nodata = dem.nodata
ew_res = 30.


# parameter value to exclude pixel
#exclude = np.nan  # i dont want use np.ma just to treat np.nan
exclude = -9999   


# factors for cardinal and diagonal directions
# for slope calculation (i.e. /length)
dx1 = 1./ew_res
dx2 = 1./(np.sqrt(2)*ew_res)

# initialize arrays
dem = np.where(dem==nodata, exclude, dem)  # exclude nodata
mask = ~np.isna(dem)
area = np.where(mask, ew_res*ns_res, exclude)
slope = np.where(mask, 0., exclude)
atb = np.where(mask, -1., exclude)


# number of pixels to calculate
natb = mask.sum() 
natbold = -1


# main loop
nmax = nrow*ncol
while((natb /= natbold) and (natb < nmax)):   
    
    # store last counter
    natbold = natb
    
    # included this print to check it it will get stuck by any reason
    print(natbold,nmax)   

    # loop through the grid cell and check if there is an upslope element
    # without an atanb value. If so, we cannot yet calculate the topidx.
    
    # ... in other words, this 1st step checks if all upslopes have been calculated.
    # which is required for the proper flow accumulation.
    for j in range(ncol):
        for i in range(nrow):      
           
            # skip non-catchment cells and cells that are done
            if( (dem[i,j] == exclude) or (atb[i,j]>0.) ):
                continue
            
            # check if river pixel
            # river cells don't accumulate flow downstream and will
            # use the average of the inflow slope as the local gradient
            if(rivermap[i,j] == 1):
                river = 1
            else:
                # not a river: check 8 flow directions for upslope elements without a topidx value ( Z[neigh]>Z[cur] & a[neigh]=none )
                not_yet = 0              
                for im in range(-1,2,1):
                    for jm range(-1,2,1):
                        ii = i+im
                        jj = j+jm
                        # check borders and ignore self pixel
                        if ( (ii >= 0 and ii < nrow) and (jj >= 0 and jj < ncol) and (im/=0 and jm/=0) ):
                            # valid pixels, check upslope and if wasnt processed
                            if ( (dem[ii,jj] /= exclude) and (dem[ii,jj] > dem[i,j]) and (atb[ii,jj]<0.) ):
                                not_yet = 1
                
                # upstream must have calculated index
                if (not_yet):
                    continue
                    
                # there are no upslope elements without a topidx value,
                # start calculations

                # find the outflow direction and calculate the sum of weights using 
                # ( tanb*countourlength). Contour length = 0.5dx for the cardinal
                # direction and 0.354dx for diagonal
                routdem = np.zeros(9) #0..8 (9 pixels)
                tanb = np.zeros(9)    #0..8 (p pixels)
                sumrt = 0
                sumtb = 0
                nrout = 0
                for im in range(-1,2,1):
                    for jm range(-1,2,1):
                        ii = i+im
                        jj = j+jm
                        if ((ii >= 0 and ii < nrow) and (jj >= 0 and jj < ncol) and (im/=0 and jm/=0) ):
                            if (dem[ii,jj] /= exclude):
                                # cardinal neighbour
                                if( (im == 0) or (jm == 0) ):
                                    dnx = dx1
                                    routefac = 0.5   #contour length
                                # diagonal neighbour
                                else:
                                    dnx = dx2
                                    routefac = 0.5 / np.sqrt(2) #contour length

                                # slope calculation when current pixel is higher than neighbour
                                if ( (dem[i,j] - dem[ii,jj]) > 0. ):
                                    tanb[k] = (dem[i,j] - dem[ii,jj]) * dnx
                                    routdem[k] = routefac * ew_res * tanb[k]  #total contour lines, weighted by (tanb*countour length).. smart.
                                    sumrt = sumrt + routdem[k]
                                    sumtb = sumtb + tanb[k]
                                    nrout = nrout + 1          #current pixel is upslope of nrout neighbours
                        # neighbour counter
                        k = k + 1
            # - end label 1
        

            # - begin sink or river cell...
            if((nrout == 0) or (river == 1)):
                nsink = nsink + 1
                river = 0

                # assume that there is a channel of length dx running midway through
                # the sink or boundary node. Take average inflow slope angle to
                # represent tanb and A/(2dx) to represent parameter 'a'
                sumtb = 0
                nslp = 0
                for im in range(-1,2,1):
                    for jm range(-1,2,1):
                        ii = i+im
                        jj = j+jm
                        if ((ii >= 0 and ii < nrow) and (jj >= 0 and jj < ncol) and (im/=0 and jm/=0) ):
                            # pathlength
                            if((im == 0) or (jm == 0)):
                                dnx = dx1
                            else:
                                dnx = dx2                            
                            if( rivermap[i,j] == 0 ):
                                # for sink/boundary squares sumtb is just the average slope (ms: for now we're summing, average comes later)
                                sumtb = subtb + (dem[ii,jj] - dem[i,j]) * dnx
                                nslp = nslp + 1
                            else:
                                # for a river cell the slope is the average inflow slope
                                if (dem[ii,jj] > dem[i,j]) :
                                    sumtb = sumtb + (dem[ii,jj] - dem[i,j]) * dnx
                                    nslp = nslp + 1
                                    
                # calculate the average inflow/sink slope angle
                if(sumtb > 0.):
                    sumtb = sumtb / nslp
                    atb[i,j] = np.log(area[i,j] / (2 * sumtb))
                    slope[i,j] = 2 * sumtb
                else:
                    atb[i,j] = exclude
                    
                # updates global counter
                natb = natb + 1
                continue                
                #- end rivers and sinks, (and current loop)


            # normal cells - previous condition for sink or "river" didnt catch 
            # average contour length
            if(sumrt > 0.):
                c = area[i,j] / sumrt    #..it look odd at first, but area is being updated in the process (see below)
                atb[i,j] = np.log(c)
                slope[i,j] = sumtb / nrout
            else: #flat neighborhood
                c = 0
                atb[i,j] = 0.01   #dummy value?
                slope[i,j] = 0
                
            # updates global counter
            natb = natb + 1
            
            # calculate downslope area: seems like an update in downslope pixels, using weighted flow accumulation
            nrout = 0
            for im in range(-1,2,1):
                for jm range(-1,2,1):
                    ii = i+im
                    jj = j+jm
                    if ((ii >= 0 and ii < nrow) and (jj >= 0 and jj < ncol) and (im/=0 and jm/=0) ):
                        if (atb[ii,jj] /= exclude):
                            if ( routdem[nrout] > 0): 
                                area[ii,jj] = area[ii,jj] + c * routdem[nrout]
                    nrout = nrout + 1 #nrout is (and must be) properly updated to match k = 0,8 in previous step
                    
# prepare output table
area = np.where(atb==exclude,exclude,area) #exclude pixels
output_atb = atb.ravel().tolist()
output_area = area.ravel().tolist()

# make dataframe, filter negative atb, save to xls
df_out = pd.DataFrame(zip(output_atb,output_area),columns=['atb','area'])
df_out_filtered = df_out[df_out['atb']>0]
df_out_filtered.to_excel('table_atb_area.xlsx')
