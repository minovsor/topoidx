'''
 Calculates the topographic index (a/ tan(beta))
  
 author: Mino Sorribas (mino.sorribas@gmail.com)
 2022.07.20
 
 Pythonized version of original C code from:
 https://github.com/ICHydro/topmodel/blob/master/topmodel/src/c_topidx.c
 
 Based on my readings I believe the C version was translated from a
 original FORTRAN code by KJ Beven & MJ Kirkby, from which i recommend
 the following paper.
  
 K. J. BEVEN & M. J. KIRKBY (1979) A physically based, variable
 contributing area model of basin hydrology / Un modèle à base physique de zone d'appel
 variable de l'hydrologie du bassin versant, Hydrological Sciences Journal, 24:1, 43-69
 DOI: 10.1080/02626667909491834

Issues:

    - I do not understand WHY rivers/sink pixels are allowed to be
    processed even if upslopes were not resolved. This is a clear
    restriction for "upland" pixels, by using the 'not_yet' variable.
    
    - River/sink pixel has the following code for the topographic index.
      atb[i,j] = np.log(area[i,j] / (2 * sumtb))
      ... but if area[i,j] must be total drainage area...
    

 
Notes:
     
    - According to Beven & Kirkby (2009) the 'a' parameter is the
     "area drained per unit contour length", but i'm not sure if people
     elsewhere are doing it properly, so here I am.
     
    - By translating the original code, it seems to me that the
    (upstream) area is updated during the process resulting in a
    flow_acc which accounts for a weighted contour length adjustment.
    (still need to run, to be sure)
    I think I'm geting close to understand what the authors had in mind.
        
    - Its not clear to me, but maybe this implementation requires a
    conditioned dem (sinks+fdr), so using raw or bare-earth dem could
    put additional challenge.  i'm too not sure, yet.
  
    - I included lots of comments to make it more understandable (to me)
    
    - Adapted running window indexes and conditions (jj,ii,etc.)

    - Replaced variable ZERO = 0.0000001  to 0. everywhere.
 
    - (Not really) issues carried from original code:
      : average downslope slope is also calculated but not returned 
      : no mechanism in place to remove negative atb values
     
 
TODO:
    - filter negative atb values and improve output table        
    - decide raster lib/reading format, probably rasterio or gdal.
    - probably include pandas for table export
    - test!
 
'''

import numpy as np


# read dem

# read river


# metadata
nrow = 
ncol = 
ew_res = 30.
nodata = -9999.

# parameters
exclude = np.nan


# factors for cardinal and diagonal directions
# for slope calculation (i.e. /length)
dx1 = 1./ew_res
dx2 = 1./(np.sqrt(2)*ew_res)

# initialize arrays
dem = np.where(dem==nodata, exclude, dem)  # exclude nodata
mask = ~np.isna(dem)
area = np.where(mask, ew_res*ns_res, np.nan)
slope = np.where(mask, 0., np.nan)
atb = np.where(mask, 0., np.nan)   # original code suggests -9.9, but why?!


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
            if( (dem[i,j] == exclude) or (atb[i,j]>0.) ):  #TODO: if stucks here, try: not np.isnan(atb[i,j])
                continue
            
            # check if river pixel
            if(rivermap[i,j] == 1):
                river = 1
            else:
                # not a river: check 8 flow directions for upslope elements without a topidx value ( Z[neigh]>Z[cur] & a[neigh]=none )
                not_yet = 0              
                for im in range(-1,2,1):
                    for jm range(-1,2,1):
                        ii = i+im
                        jj = j+jm
                        if ( (ii >= 0 and ii < nrow) and (jj >= 0 and jj < ncol) and (im/=0 and jm/=0) ): # respect matrix borders and ignore self pixel
                            if ( (dem[ii,jj] /= exclude) and (dem[ii,jj] > dem[i,j]) and (np.isnan(atb[ii,jj])) ):  #
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
                c = area[i,j] / sumrt    #..it look odd at first, but area is  being updated in the process (see below)
                atb[i,j] = np.log(c)
                slope[i,j] = sumtb / nrout
            else:
                c = 0
                atb[i,j] = 0.01
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
                    
  # format output
  output_atb = []
  output_area = []
  for j in range(ncol):
    for i in range(nrow):
        if ( atb[i,j] == exclude):
            area[i,j] = exclude      
        output_atb = atb[i,j]
        #output_area = area[i,j] #not rly needed
