######Data Module######
#VJS 6/2016

#Module to read and digest data of different forms, and prepare it for grmpepy

def mread(flatfile,hashfile,stationfile,station_cols):
    '''
    Read data from Annemarie's flatfile
    VJS 6/2016
    
    Input 
        flatfile:        String with path to the anza flatfile from AB
        hashfile:        String with path to the anza hash file from AB, with locations
        stationfile:     STring with path to the station file, for anza stations 
        station_cols: 	 Array with columns to use for station name, lat, lon, and el [name, lat, lon, el]
    Output
        event:          Array with event numbers
        sta:            List with station names
        N:              Array with station numbers (stnum)
        ml:             Array with local magnitudes
        mw:             Array with moment magnitudes
        DA:             Array with PGA values, geometrically averaged from two components, in nm/s/s
        DV:             Array with PGV values, geometrically averaged from two components, in nm/s/s
        dist:           Array with epicentral distance, km
        vs30:           Array with vs30 values
        lat:            Array with event latitude
        lon:            Array with event longitude
        depth:          Array with event depth
        stlat:          Array with station latitude
        stlon:          Array with station longitude
        stelv:          Array with station elevation
        source_i:       Array with the source index for raytracing (i.e., unique sources numbered)
        receiver_i:     Array with the receiver index for raytracing (i.e., unique receivers numbered)
    '''
    
    import scipy.io as sio
    import cdefs
    from numpy import genfromtxt,where,zeros,unique
    
    
    #Read in flatfile
    datr=sio.loadmat(flatfile)
    hashr=genfromtxt(hashfile)
    #Read in station info from columns of station file:
    #Station name:
    name_col=station_cols[0]
    lat_col=station_cols[1]
    lon_col=station_cols[2]
    elv_col=station_cols[3]
    stat_sta_r=genfromtxt(stationfile,skip_header=3,usecols=name_col,dtype='S5')
    stat_coord_r=genfromtxt(stationfile,skip_header=3,usecols=[lat_col,lon_col,elv_col])
    
    #Info from the hashfile:
    hevent=hashr[:,6]
    hlat=hashr[:,7]
    hlon=hashr[:,8]
    hdepth=hashr[:,9]
    
    #Extract components:
    devent=datr['event']            #event number
    dsta=datr['sta']                #station name
    dhdrs=datr['hdrs']              #header info...see below...
    dN=datr['N']                    #station number         
    dMl=datr['Ml']                  #catalog/local mag.
    dMw=datr['Mw']                  #Ml converted to Mw
    dDA=datr['DA']                  #PGA, geometrically averaged from 2 horiz. comp.
    dDV=datr['DV']                  #PGV, geometrically averaged from 2 horiz. comp.
    dPGA=datr['PGA']                #DA, adjusted for Vs30 and converted to g
    dPGV=datr['PGV']                #DV, adjusted for Vs30 and converted to g
    dlogPGA=datr['logPGA']          #DA adjusted to g
    dlogPGV=datr['logPGV']          #DV adjusted to cm/s
    dnewlogPGA=datr['newlogPGA']    #log10 of PGA adjusted to 10km
    dnewlogPGV=datr['newlogPGV']    #log10 of PGV adjusted to 10km
    dVs30=datr['Vs30']              #Vs30 for each station, reference in paper.
    
    #Find the lat/lon that corresponds to this event:
    #Zero out lat/lon arrays:
    ev_lat=zeros((len(devent[0]),1))
    ev_lon=zeros((len(devent[0]),1))
    ev_dep=zeros((len(devent[0]),1))
    
    ######
    #Get data from hashfile, put it in here...
    #Find which row (called event_ind) in the hashfile corresponds to this event:
    for i in range(len(devent[0])):
        event_ind=where(hevent==devent[0,i])[0][0]
        #Set that row in lat/lon to this value:
        ev_lat[i]=hlat[event_ind]
        ev_lon[i]=hlon[event_ind]
        ev_dep[i]=hdepth[event_ind]
    
    #Find row from station data that corresponds to the station:
    #First zero out arrays:
    sta_lat=zeros((len(dsta[0]),1))
    sta_lon=zeros((len(dsta[0]),1))
    sta_elv=zeros((len(dsta[0]),1))
    
    for j in range(len(dsta[0])):
        stationi=dsta[0][j][0].astype('S5')
        #Find which index of the station file this corresponds to for this station:
        station_ind=where(stationi==stat_sta_r)[0]
        
        #Now get the information for this station, lat and lon:
        sta_lat[j]=stat_coord_r[station_ind,0]
        sta_lon[j]=stat_coord_r[station_ind,1]
        sta_elv[j]=stat_coord_r[station_ind,2]
        
    
    ###Get indices for event and station for the sources.in and receivers.in files####
    ##Events first:
    
    #Get the unique station and event indices:
    unique_events=unique(devent)
    
    #Zero out source ind array:
    source_ind=zeros((len(devent[0])))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(devent[0]==eventi)
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
    
    #Now set these to integers...
    source_ind=source_ind.astype('int64')
    
    ##Next stations:
    unique_stations=unique(dsta)
    
    #Zero out array:
    receiver_ind=zeros((len(dsta[0])))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(dsta[0]==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_ind=receiver_ind.astype('int64')
  
        
    #Put into arrays (instead of array of arrays):
    event=devent[0]
    sta=dsta[0]
    N=dN[0]
    ml=dMl[:,0]
    mw=dMw[:,0]
    DA=dDA[:,0]
    DV=dDV[:,0]
    PGA=dPGA[:,0]
    PGV=dPGV[:,0]
    logPGA=dlogPGA[:,0]
    logPGV=dlogPGV[:,0]
    newlogPGA=dnewlogPGA[:,0]
    newlogPGV=dnewlogPGV[:,0]
    vs30=dVs30[0]
    lat=ev_lat[:,0]
    lon=ev_lon[:,0]
    depth=ev_dep[:,0]
    stlat=sta_lat[:,0]
    stlon=sta_lon[:,0]
    stelv=sta_elv[:,0]
    source_i=source_ind
    receiver_i=receiver_ind
    
    ##
    #Split apart the header (dhdrs)
    ddepmin=dhdrs['depmin']         #dependent (y) min
    ddepmax=dhdrs['depmax']         #dependent (y) max
    dmag=dhdrs['mag']               #magnitude (catalog)
    ddist=dhdrs['dist']             #distance from source to site, epicentral...
    daz=dhdrs['az']                 #source-site azimuth
    dyear=dhdrs['year']             #year
    dday=dhdrs['day']               #day
    dhour=dhdrs['hour']             #hour
    dmin=dhdrs['min']               #min
    dsec=dhdrs['sec']               #sec
    dmsec=dhdrs['msec']             #msec
    dnevid=dhdrs['nevid']           #event number
    didep=dhdrs['idep']             #units of the independent variable:
                                    #5 implies no instrument response removed
                                    #7 is velocity in nm/s
                                    #8 is acceleration in nm/s/s
    
    #Put into useful values...
    depmin=ddepmin[0][0]
    depmax=ddepmax[0][0]
    mag=dmag[0][0]
    dist=ddist[0][0]
    az=daz[0][0]
    year=dyear[0][0]
    day=dday[0][0]
    minu=dmin[0][0]
    sec=dsec[0][0]
    msec=dmsec[0][0]
    nevid=dnevid[0][0]
    idep=didep[0][0]
    
    #Return the event numeber, station name, station number, ml,mw, PGA,pgv, 
    #epcentral distance (Dist), vs30
    return event,sta,N,ml,mw,DA,DV,dist[:,0],vs30,lat,lon,depth,stlat,stlon,stelv,source_i,receiver_i



#############################################################################
def find_source_receiver_indices(evnum,sta):
    '''
    Given a set of events and stations from a database of recordings, 
    not ordered in their arrays to be unique, what are their unique
    indices in the same given order?
    Input:
        evnum:          Array from database of recordings with event number for each recording
        sta:            Array from database of recordings with station name for each recording
    Output:
        source_i:       Array with the unique identifier (starting at 1) for the sources for each recording
        receiver_i:     Array with the unique identifier (starting at 1) for the receivers for each recording
    '''
    
    from numpy import unique,zeros,where
    
    ###Get indices for event and station for the sources.in and receivers.in files####
    ##Events first:
    
    #Get the unique station and event indices:
    unique_events=unique(evnum)
    
    #Zero out source ind array:
    source_ind=zeros((len(evnum)))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(evnum==eventi)
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
    
    #Now set these to integers...
    source_ind=source_ind.astype('int64')
    
    ##Next stations:
    unique_stations=unique(sta)
    
    #Zero out array:
    receiver_ind=zeros((len(sta)))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(sta==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_ind=receiver_ind.astype('int64')
    
    ######
    #At the end, convert the lists to arrays:
    source_i=source_ind
    receiver_i=receiver_ind
    
    return source_i, receiver_i


#############################################################################

def read_obj_list(objfile):
    '''
    Read a pickle file with a list of event of station objects
    Input:
        objfile:        String with path to the pickle file containing the list of objects
    Output:
        obj_list:   List with event or station objects in it
    '''
    
    import cPickle as pickle
    
    #Open the file
    obj=open(objfile,'r')
    
    #Zero out a list to add the event objects to:
    obj_list=[]
    
    #REad the data...
    readfile=True
    while readfile==True:
        try:
            obji=pickle.load(obj)
            obj_list.append(obji)
        except EOFError:
            print ('File read complete - '+objfile)
            readfile=False
            obj.close()
        
    return obj_list
    
    
    
######
def db_station_sample(dbpath_in,numstas,dbpath_out):
    '''
    Sample a database to only include events recorded on a minimum number
    of stations
    VJS 8/2016
    
    Input:
        dbpath_in:          String with path to the input database
        numstas:            Minimum number of stations for events to be recorded on
        dbpath_out:         STring with path to output database
    Output:
        Writes out the sampled database to dbpath_out
    '''
    import cPickle as pickle
    from numpy import unique,where,array,r_,zeros
    import cdefs as cdf
    
    #Read in the original database
    dbfile=open(dbpath_in,'r')
    db_orig=pickle.load(dbfile)
    dbfile.close()
    
    #Find how many unique events there are:
    unique_events=unique(db_orig.evnum)
    nevents=len(unique_events)
    
    #Initiate the "keep" event index array:
    keep_event_ind=array([]).astype('int')
    
    #Loop over the unique events and determine if they are recorded on the 
    #minimum number of stations:
    
    for event_ind in range(nevents):
        #Call the event:
        eventi=unique_events[event_ind]
        #Find where in the database there are recordings of this event:
        db_event_ind=where(db_orig.evnum==eventi)[0].astype('int')
        
        #Get the stations for this event:
        eventi_stas=db_orig.sta[db_event_ind]
        #Get the unique stations recording this event:
        num_unique_stas_i=len(eventi_stas)
        
        #If it's at least the number of minimum stations, keep this stuff, save 
        #it to the keep index variable:
        if num_unique_stas_i>=numstas:
            keep_event_ind=r_[keep_event_ind,db_event_ind]
    
    
    #Save the keep event variable as integers:
    keep_event_ind.astype('int')
    
    #Now save just these indices in the database:
    edepth=db_orig.edepth[keep_event_ind]
    elat=db_orig.elat[keep_event_ind]
    elon=db_orig.elon[keep_event_ind]
    evnum=db_orig.evnum[keep_event_ind]
    ffdf=db_orig.ffdf[keep_event_ind]
    md_ffdf=db_orig.md_ffdf[keep_event_ind]
    ml=db_orig.ml[keep_event_ind]
    mw=db_orig.mw[keep_event_ind]
    pga=db_orig.pga[keep_event_ind]
    pga_pg=db_orig.pga_pg[keep_event_ind]
    pgv=db_orig.pgv[keep_event_ind]
    r=db_orig.r[keep_event_ind]
    sta=db_orig.sta[keep_event_ind]
    stlat=db_orig.stlat[keep_event_ind]
    stlon=db_orig.stlon[keep_event_ind]
    stelv=db_orig.stelv[keep_event_ind]
    stnum=db_orig.stnum[keep_event_ind]
    vs30=db_orig.vs30[keep_event_ind]
    
    if db_orig.vs30_method!=None:
        vs30_method=db_orig.vs30_method[keep_event_ind]
    else:
        vs30_method=None
        
    if db_orig.pga_snr!=None:
        pga_snr = db_orig.pga_snr[keep_event_ind]
    else:
        pga_snr = None
        
    if db_orig.pgv_snr!=None:
        pgv_snr = db_orig.pgv_snr[keep_event_ind]
    else:
        pgv_snr = None
    
    ###Change the source and receiver indices for raytracing...
    #Get the unique station and event indices:
    unique_events=unique(evnum)
    
    #Zero out source ind array:
    source_ind=zeros((len(evnum)))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(evnum==eventi)
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
        
    #Now set these to integers...
    source_i=source_ind.astype('int64')
    
    ##
    ##Next stations:
    unique_stations=unique(stnum)
    
    #Zero out array:
    receiver_ind=zeros((len(stnum)))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(stnum==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_i=receiver_ind.astype('int64')
    
    ##BEFORE SAVING:
    ##cdefs only takes DA and DV in nm/s/s and nm/s...convert to these (currently
    ##in m/s/s and m/s)
    #DA=pga/1e-9
    #DV=pga/1e-9
    
    #Make sampled database:
    db_samp=cdf.db(evnum,sta,stnum,ml,mw,pga,pgv,r,vs30,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i,vs30_method=vs30_method,pga_snr=pga_snr,pgv_snr=pgv_snr)
    
    #Save to file...
    doutfile=open(dbpath_out,'w')
    pickle.dump(db_samp,doutfile)
    doutfile.close()
    
    
######
def recording_sample(dbpath_in,recording_indices,dbpath_out):
    '''
    Sample a database to only include certain indices
    of stations
    VJS 8/2016
    
    Input:
        dbpath_in:          String with path to the input database
        recording_indices:  Array with the indices of recordings to KEEP
        dbpath_out:         STring with path to output database
    Output:
        Writes out the sampled database to dbpath_out
    '''
    import cPickle as pickle
    from numpy import unique,where,array,r_,zeros,setdiff1d
    import cdefs as cdf
    
    #Read in the original database
    dbfile=open(dbpath_in,'r')
    db_orig=pickle.load(dbfile)
    dbfile.close()
    
    #Find how many unique events there are:
    unique_events_orig=unique(db_orig.evnum)
    nevents_orig=len(unique_events_orig)
    
    unique_sta_orig=unique(db_orig.sta)
    nsta_orig=len(unique_sta_orig)
    
    #Now save just these indices of the "keep" recordings in the database:
    edepth=db_orig.edepth[recording_indices]
    elat=db_orig.elat[recording_indices]
    elon=db_orig.elon[recording_indices]
    evnum=db_orig.evnum[recording_indices]
    ffdf=db_orig.ffdf[recording_indices]
    md_ffdf=db_orig.md_ffdf[recording_indices]
    ml=db_orig.ml[recording_indices]
    mw=db_orig.mw[recording_indices]
    pga=db_orig.pga[recording_indices]
    pga_pg=db_orig.pga_pg[recording_indices]
    pgv=db_orig.pgv[recording_indices]
    r=db_orig.r[recording_indices]
    sta=db_orig.sta[recording_indices]
    stlat=db_orig.stlat[recording_indices]
    stlon=db_orig.stlon[recording_indices]
    stelv=db_orig.stelv[recording_indices]
    stnum=db_orig.stnum[recording_indices]
    vs30=db_orig.vs30[recording_indices]
    
    if db_orig.vs30_method!=None:
        vs30_method=db_orig.vs30_method[recording_indices]
    else:
        vs30_method=None
        
    if db_orig.pga_snr!=None:
        pga_snr = db_orig.pga_snr[recording_indices]
    else:
        pga_snr = None
        
    if db_orig.pgv_snr!=None:
        pgv_snr = db_orig.pgv_snr[recording_indices]
    else:
        pgv_snr = None
    
    ###Change the source and receiver indices for raytracing...
    #Get the unique station and event indices:
    unique_events=unique(evnum)
    
    #Zero out source ind array:
    source_ind=zeros((len(evnum)))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(evnum==eventi)
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
        
    #Now set these to integers...
    source_i=source_ind.astype('int64')
    
    ##
    ##Next stations:
    unique_stations=unique(stnum)
    
    #Zero out array:
    receiver_ind=zeros((len(stnum)))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(stnum==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_i=receiver_ind.astype('int64')
    
    ##BEFORE SAVING:
    ##cdefs only takes DA and DV in nm/s/s and nm/s...convert to these (currently
    ##in m/s/s and m/s)
    #DA=pga/1e-9
    #DV=pga/1e-9
    
    #Make sampled database:
    db_samp=cdf.db(evnum,sta,stnum,ml,mw,pga,pgv,r,vs30,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i,vs30_method=vs30_method,pga_snr=pga_snr,pgv_snr=pgv_snr)
    
    #Save to file...
    doutfile=open(dbpath_out,'w')
    pickle.dump(db_samp,doutfile)
    doutfile.close()
    
    # Print stats:
    print('Originally %i events, now %i events' % (nevents_orig,len(unique(evnum))))
    print('Originally %i stations, now %i stations' % (nsta_orig,len(unique(sta))))
    
    #   Difference of old and new event sets:
    event_diff=setdiff1d(unique_events_orig,unique(evnum))
    sta_diff=setdiff1d(unique_sta_orig,unique(sta))
    
    print('The events removed are: \n')
    print(event_diff)
    print('\n The stations removed are: \n')
    print(sta_diff)
    
    
def db_propgrid_sample(dbpath_in,propgrid,dbpath_out):
    '''
    Sample a database to only include events recorded on a minimum number
    of stations
    VJS 8/2016
    
    Input:
        dbpath_in:          String with path to the input database
        propgrid:           Propagation grid limits in format: [[W,E],[S,N]]
        dbpath_out:         STring with path to output database
    Output:
        Writes out the sampled database to dbpath_out
    '''
    import cPickle as pickle
    from numpy import unique,where,array,r_,zeros
    import cdefs as cdf
    from matplotlib import path
    
    #Read in the original database
    dbfile=open(dbpath_in,'r')
    dbin=pickle.load(dbfile)
    dbfile.close()
    
    # Make the propgrid outline path:
    w = propgrid[0][0]
    e = propgrid[0][1]
    s = propgrid[1][0]
    n = propgrid[1][1]
    
    gridpath = path.Path([[w,s],[e,s],[e,n],[w,n],[w,s]])
    
    print('Read file %s' % dbpath_in)
    print('%s recordings read in, with %s unique events, and %s unique stations' % (str(len(dbin.evnum)),str(len(unique(dbin.evnum))),str(len(unique(dbin.sta)))))
    
    
    # Loop through the recordings, and if both the event and station are inside
    #   the path, keep the recording:
    
    #Initiate the "keep" event index array:
    keep_event_ind=array([]).astype('int')
    
    for record_i in range(len(dbin.evnum)):
        if (gridpath.contains_point([dbin.elon[record_i],dbin.elat[record_i]])) & (gridpath.contains_point([dbin.stlon[record_i],dbin.stlat[record_i]])):
           
            keep_event_ind = r_[keep_event_ind,record_i]
            
    evnum=dbin.evnum[keep_event_ind]
    sta=dbin.sta[keep_event_ind]
    stnum=dbin.stnum[keep_event_ind]
    ml=dbin.ml[keep_event_ind]
    mw=dbin.mw[keep_event_ind]
    pga=dbin.pga[keep_event_ind]
    pgv=dbin.pgv[keep_event_ind]
    pga_pg=dbin.pga_pg[keep_event_ind]
    r=dbin.r[keep_event_ind]
    vs30=dbin.vs30[keep_event_ind]
    ffdf=dbin.ffdf[keep_event_ind]
    md_ffdf=dbin.md_ffdf[keep_event_ind]
    elat=dbin.elat[keep_event_ind]
    elon=dbin.elon[keep_event_ind]
    edepth=dbin.edepth[keep_event_ind]
    stlat=dbin.stlat[keep_event_ind]
    stlon=dbin.stlon[keep_event_ind]
    stelv=dbin.stelv[keep_event_ind]
    
    if dbin.vs30_method!=None:
        vs30_method=dbin.vs30_method[keep_event_ind]
    else:
        vs30_method=None
    
    if dbin.pga_snr!=None:
        pga_snr = dbin.pga_snr[keep_event_ind]
    else:
        pga_snr = None
        
    if dbin.pgv_snr!=None:
        pgv_snr = dbin.pgv_snr[keep_event_ind]
    else:
        pgv_snr = None
    
    
    
    print('Data reduced to %s unique events, and %s unique stations' % (str(len(unique(evnum))),str(len(unique(sta)))))
        
        
        
    # Now get source_i and receiver_i for the new dataset:
    
    ###Change the source and receiver indices for raytracing...
    #Get the unique station and event indices:
    unique_events=unique(evnum)
    
    #Zero out source ind array:
    source_ind=zeros((len(evnum)))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(evnum==eventi)[0]
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
        
    #Now set these to integers...
    source_i=source_ind.astype('int64')
    
    ##
    ##Next stations:
    unique_stations=unique(stnum)
    
    #Zero out array:
    receiver_ind=zeros((len(stnum)))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(stnum==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_i=receiver_ind.astype('int64')
    
    
    ## Now make new sampled database:
    dbsamp = cdf.db(evnum,sta,stnum,ml,mw,pga,pgv,r,vs30,elat,elon,edepth,stlat,stlon,stelv,source_i,receiver_i,vs30_method=vs30_method,pga_snr=pga_snr,pgv_snr=pgv_snr)
    
    #Save to file...
    doutfile=open(dbpath_out,'w')
    pickle.dump(dbsamp,doutfile)
    doutfile.close()
    
    
def multiseg2pckl(multisegpath,pcklpath,pathlimits):
    '''
    VJS 9/2016
    Convert a GMT multisegment file to a pckl file to be plotted in python
    Input:
        multisegpath:       String with the path to the input multisegment file
        pcklpath:           String with the path to the output pckl file; List 
                            of arrays, each with a segment to scatter or plot
                            Output to pcklpath.
        pathlimits:         [[lonmin,lonmax],[latmin,latmax], same as
                            [[xmin,xmax],[ymin,ymax]] 
    Output:
        allsegments         List of arrays, each with two columns (lon, lat)
    '''
    
    from numpy import zeros,array,where
    import matplotlib.path as mplPath
    import cPickle as pickle
    
    #Get the corner coordinates out of pathlimits: first lonmin/max, latmin/max:
    lonmin=pathlimits[0][0]
    lonmax=pathlimits[0][1]
    latmin=pathlimits[1][0]
    latmax=pathlimits[1][1]
    
    #Now the actual corners:
    bottom_left=[lonmin,latmin]
    bottom_right=[lonmax,latmin]
    top_right=[lonmax,latmax]
    top_left=[lonmin,latmax]
    
    #Define the path bounds - for a squre, it's 5 points:
    path_coordinates=array([bottom_left,bottom_right,top_right,top_left,bottom_left])
    #Define the regionpath with the mplPath command:
    region_path=mplPath.Path(path_coordinates)
    
    #First count the number of segments and the number of elements in each segment
    Nsegments=0
    Num_elements=[]
    f=open(multisegpath,'r')
    first_line=True
    while True:
        line=f.readline()
        if '>' in line:
            if first_line==False: #Append previous count
                Num_elements.append(Numel)
            first_line=False
            Nsegments+=1
            Numel=0
        else:
            Numel+=1
        if line=='': #End of file
            Num_elements.append(Numel-1)
            break
    f.close()
            
    #Now loop over segments and make an arra per segment adn append to list of arrays
    all_segments=[]    
    f=open(multisegpath,'r')
    
    for ksegment in range(Nsegments):
        
        #First line in the segmetn is the stupid carrot
        line=f.readline()
        
        #Now read the next Num_elements[ksegment] lines
        lonlat=zeros((Num_elements[ksegment],2))
        for kelement in range(Num_elements[ksegment]):
            line=f.readline()
            lonlat[kelement,0]=float(line.split()[0])
            lonlat[kelement,1]=float(line.split()[1])
            
        
        #Before appending this segment to the list, check if any points along 
        #the segment are in the path
        
        #Are any points of this segment in the path defined above?
        points_logical=where(region_path.contains_points(lonlat)==True)[0]
        
        #If any points along htis segment are contained:
        if len(points_logical>0):  
            #Done, append to list
            all_segments.append(lonlat)
    
    f.close()
    
    #Write to the pickle file:
    fout=open(pcklpath,'w')
    for segment_i in range(len(all_segments)):
        pickle.dump(all_segments[segment_i],fout)
    fout.close()
    
    return all_segments
                
    
    
def read_material_model(coordspath,modelpath):
    '''
    Read in a velocity model (like Fang 2016), parse into format to be read by
    cdefs to make an object
    Input:
        coordspath:             String with path to the coordinates file, with info
                                    about the x, y, and z limits of the model
                                    File format:
                                        x1 x2 x3 ...\r\n
                                        \r\n
                                        y1 y2 y3 ...\r\n
                                        \r\n
                                        z1  z3  z3  ...   (double spaces between z's)
                                        
        modelpath:           String with path to the velocity file (i.e., Vp or Vs)
                                    with format columns: x, rows: y; repeats in z
    Output:
        x:                      Array with the x values of the model nodes
        y:                      Array with the y values of the model nodes
        z:                      Array with the z values of the model nodes
        nx:                     Float with number of x points
        ny:                     Float with number of y points
        nz:                     Float with number of z points
        model:                  Multi-dim array with model: len(model) = len(z);
                                    shape(model[0]) = len(y),len(x)
    '''
    from numpy import array,zeros,genfromtxt,shape
    
    #Read in the model info and store:
    model_in=genfromtxt(modelpath)
    
    ###Read in the coordinates file:
    cfile=open(coordspath,'r')
    
    #Read line by line, starting with x:
    xline_raw=cfile.readline()
    #split up - first cut off the \r\n, then split each number by spaces.
    xline=array(xline_raw.split('\r')[0].split(' '))
    x=xline.astype(float)
    
    #Read one more line since there's a space:
    cfile.readline()
    
    #Again for y and z:
    yline_raw=cfile.readline()
    #split up - first cut off the \r\n, then split each number by spaces.
    yline=array(yline_raw.split('\r')[0].split(' '))
    y=yline.astype(float)
     
    cfile.readline()                   
                    
    zline_raw=cfile.readline()
    #split up - first cut off the \r\n, then split each number by spaces.
    zline=array(zline_raw.split('  '))
    z=zline.astype(float)         
       
    #Close file:
    cfile.close()   
    
                    
    ###########
    ##Now get model info:
    #Number of x, y, and z points:
    nx=len(x)
    ny=len(y)
    nz=len(z)
    
    #Initialize the model array:
    material_model=zeros((nz,ny,nx))
                
    #Now extract model values.
    #Loop over the number of z entries, and pull out the chunk that corresponds 
    #to that z; then add it to material_model
    #Initiate a counter for counting the z position:
    count_z=0
    
    for z_i in range(nz):
        i_beg=z_i
        i_end=z_i+ny
        material_model[z_i]=model_in[count_z:count_z+ny,:]
        
        #Add to counter:
        count_z=count_z+ny
            
    
    #Return values:
    return x, y, z, nx, ny, nz, material_model


############################################################################
def read_hauksson_file(Qtxtpath,provide_simple=False):
    '''
    Read in Egill's Q model (Hauksson and Shearer 2006), parse into format to be read by
    cdefs to make an object
    Input:
        Qtxtpath:             String with path to the model file, with info
                                    about the x, y, and z limits of the model
                                    and Q data
                                    File format:
                                        LAYER 2             1.00km      5.36kmsec-1
                                         long., lat., percent velo change, abs. velocity
                                         QpLatdegrLatminu  LondegrLonminu  -X-Ygrid  Depth(km)
        provide_simple:       String with flag to provide a simple output of only the points of the model, in form: modelx,modely,modelz,modeldata
                                    where each are vectors of length n, n being the total number of nodes in the model
    Output:
        lon:                      Array with the x values of the model nodes
        lat:                      Array with the y values of the model nodes
        depth:                    Array with the z values of the model nodes
        nx:                       Float with number of unique x (longitude) values 
        ny:                       Float with number of unique y (latutide) values
        nz:                       Float with number of unique z (depth) values
        Q:                        Multi-dim array with model: len(model) = len(z);
                                    shape(model[0]) = len(y),len(x)
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    
    # Open file:
    f = open(Qtxtpath,'r')

    x = []
    y = []
    z = []
    
    # Now with the number of layers:
    f = open(Qtxtpath,'r')
    layercounter = 0
    layer_lines = []

    print('Reading in lines to get the number of layers')
    for line in f.readlines():
        if 'LAYER' in line:
            layercounter+=1
            print(line)
            if layercounter == 1:
                i_layer_line = []
            else:
                print('appending')
                layer_lines.append(i_layer_line)
                i_layer_line = []
        if ('LAYER' not in line) and ('long.' not in line):
                i_layer_line.append(line)
        
    # Append last one:
    layer_lines.append(i_layer_line)
    
    f.close()
            
    # Now that all are in the list, loop through and pull out the info for each
    #  z-slice:
    
    print(' Initializing arrays to parse information...')
    
    lon_deg = []
    lon_min = []
    lat_deg = []
    lat_min = []
    depth_arr = []
    Q_arr = []

    # Final lon and lat arrays:
    lon_arr = []
    lat_arr = []
    
    print('Starting to loop over layers to get information')
    # Loop over the layers:
    for i_layer in range(len(layer_lines)):
        ilayer_lon_deg = []
        ilayer_lon_min = []
        ilayer_lat_deg = []
        ilayer_lat_min = []
        ilayer_depth_arr = []
        ilayer_Q_arr = []

        # Loop over the strings in this layer:
        i_layer_line = layer_lines[i_layer]

        for j_line in range(len(i_layer_line)):
            ij_line = i_layer_line[j_line]
            
            ij_Q = np.float(ij_line[3:10])
    
            ij_lat_deg = np.float(ij_line[10:12])
            ij_lat_min = np.float(ij_line[12:16])
    
            ij_lon_deg = np.float(ij_line[18:21])
            ij_lon_min = np.float(ij_line[21:26])
    
            ij_depth = np.float(ij_line[52:57])
            
            # Append these to this layer's lists:
            ilayer_lat_deg.append(ij_lat_deg)
            ilayer_lat_min.append(ij_lat_min)
            ilayer_lon_deg.append(ij_lon_deg)
            ilayer_lon_min.append(ij_lon_min)
            ilayer_depth_arr.append(ij_depth)
            ilayer_Q_arr.append(ij_Q)
        
        # Turn these layer lists into arrays:
        ilayer_lat_deg = np.array(ilayer_lat_deg)
        ilayer_lat_min = np.array(ilayer_lat_min)
        ilayer_lon_deg = np.array(ilayer_lon_deg)
        ilayer_lon_min = np.array(ilayer_lon_min)
        ilayer_depth_arr = np.array(ilayer_depth_arr)
        ilayer_Q_arr = np.array(ilayer_Q_arr)

        # Append this layer's arrays to the master lists:       
        lat_deg.append(ilayer_lat_deg)
        lat_min.append(ilayer_lat_min)
        lon_deg.append(ilayer_lon_deg)
        lon_min.append(ilayer_lon_min)
        depth_arr.append(ilayer_depth_arr)
        Q_arr.append(ilayer_Q_arr)
        
        # Convert lon/lat to decimal degrees:
        ilayer_lat_arr = np.round(ilayer_lat_deg + (ilayer_lat_min/60.),decimals=6)
        ilayer_lon_arr = -1.*np.round(ilayer_lon_deg + (ilayer_lon_min/60.),decimals=6)

        # Append to main lon/lat arrays:
        lat_arr.append(ilayer_lat_arr)
        lon_arr.append(ilayer_lon_arr)
        
    # Now, loop through layers.  For each layer, findthe number of unique lats
    #    and lons.  It should be the same for each depth slice:
    i_num_uniquelon_list = []
    i_num_uniquelat_list = []
    for ilayer in range(len(layer_lines)):
        i_num_uniquelon = len(np.unique(lon_arr[ilayer]))
        i_num_uniquelon_list.append(i_num_uniquelon)
        
        i_num_uniquelat = len(np.unique(lat_arr[ilayer]))
        i_num_uniquelat_list.append(i_num_uniquelat)
        
    # At end find number that are the same:
    uniquelon_count = len(np.unique(i_num_uniquelon))
    uniquelat_count = len(np.unique(i_num_uniquelat))
    
    if (uniquelon_count == 1) & (uniquelat_count == 1):
        print('Same number of unique lon and lat points in each layer')
        nx = i_num_uniquelon_list[0]
        ny = i_num_uniquelat_list[0]
        nz = len(layer_lines)
        
        unique_x = np.unique(lon_arr[0])
        unique_y = np.unique(lat_arr[0])
        unique_z = np.unique(depth_arr)
        
    ########
    # Get grid for data:
    gridx = np.sort(unique_x)
    gridy = np.sort(unique_y)
    gridX,gridY = np.meshgrid(gridx,gridy)
        
    Q_grid = np.zeros((nz,ny,nx))
    
    for i_depth in range(nz):
        i_existing_points = np.c_[lon_arr[i_depth],lat_arr[i_depth]]
        i_Q_slice = griddata(i_existing_points,Q_arr[i_depth],(gridX,gridY),method='cubic',fill_value=0)
        Q_grid[i_depth] = i_Q_slice
        
    ## If provide_simple was set to False, only return this infor:
    if provide_simple == False:
        return unique_x, unique_y, unique_z, nx, ny, nz, Q_grid
        
    ## OTherwise, if provide_simple is true, then provide lon_arr, lat_arr, depht_arr, and Q_arr in raveled formats, just arrays:
    elif provide_simple == True:
        nodex = np.array(lon_arr).ravel()
        nodey = np.array(lat_arr).ravel()
        nodez = np.array(depth_arr).ravel()
        nodeQ = np.array(Q_arr).ravel()
        
        # Return values:
        return unique_x, unique_y, unique_z, nx, ny, nz, Q_grid, nodex, nodey, nodez, nodeQ
        

#####
#Read in Janine's PGA format file
def read_jsbfile(datafile,get_year='no'):
    '''
    Read in Janine's PGA/PGV data format file and print out usable data for the db object read.
    Input:
        datafile:           String with path to the datafile
        get_year:           String to include/output year - 'yes' or 'no'
    Output:
        evnum:              Array with event numbers
        evlat:              Array with event latitude
        evlon:              Array with event longitude
        evdep:              Array with event depth (depth positive)
        sta:                Array with station name
        stlat:              Array with station latitude
        stlon:              Array with station longitude
        stelv:              Array with station elevation (elevation positive)
        grcircle:           Array with great circle paths
        ml:                 Array with local magnitudes
        mw:                 Array with moment maginitudes
        pga_mgal:           Array with PGA in milligals
        source_i:           Array with the source number for each recoridng, for raytracing
        receiver_i:         Array with the receiver number for each recoridng, for raytracing
        predparam_snr:      Array with signal to noise ratio for every recording
        evyear:             Array with event year, if requested
    '''
    
    from numpy import genfromtxt,unique,log10,array,where,zeros
    
    #First read in the stations from the flatfile:
    #Will be in two columns:   sta    chan
    sta_r=genfromtxt(datafile,dtype='S5',usecols=[0])
    chan_r=genfromtxt(datafile,dtype='S5',usecols=[1])

    #Read in the data from the flatfile:
    dat_r=genfromtxt(datafile,usecols=range(2,14))
    
    #Set variables from dat_r:
    stlat_r=dat_r[:,0]
    stlon_r=dat_r[:,1]
    stelv_r=dat_r[:,2]
    evnum_r=dat_r[:,3]
    evlat_r=dat_r[:,4]
    evlon_r=dat_r[:,5]
    evdep_r=dat_r[:,6]
    grcircle_r=dat_r[:,7]
    ml_r=dat_r[:,8]
    mw_r=dat_r[:,9]
    # If this is pga, units are mgal; if it's pgv, units are cm/s
    predparam_r=dat_r[:,10]
    predparam_snr_r=dat_r[:,11]
    
    # If year is requested, grab event origin date:
    if get_year=='yes':
        evdate_r = genfromtxt(datafile,dtype='S',usecols=[15])
    
    #Set the final variables to empty lists, to append to:
    evnum=[]
    evlat=[]
    evlon=[]
    evdep=[]
    sta=[]
    stlat=[]
    stlon=[]
    stelv=[]
    grcircle=[]
    ml=[]
    mw=[]
    # If it's PGA this is mgal; if it's PGV this is cm/s.
    predparam=[]
    predparam_snr=[]
    
    if get_year=='yes':
        evdate=[]
    
    #Find the geometrical average for the E and N stations of each event:
    #First get the unique events:
    unique_events=unique(evnum_r)
    
    #For every unique event, find which stations record it:
    for event_i in range(len(unique_events)):
        unique_event_ind=where(evnum_r==unique_events[event_i])[0]
        #Stations recoridng this event:
        sta_event_iter=sta_r[unique_event_ind]
        #Get the unique stations:
        unique_stations_iter=unique(sta_event_iter)
        
        #For each station, get the E and N component and average them:
        for station_i in range(len(unique_stations_iter)):
            #What are the indices of this station in sta_r, same length as chan_r,
            #for both E and N?
            unique_station_ind=where(sta_r[unique_event_ind]==unique_stations_iter[station_i])[0]
            
            #Set the channel E and N counter:
            E_counter=0
            N_counter=0
            #Loop through the channels and get E and N:
            for chan_iter in range(len(unique_station_ind)):
                channel=chan_r[unique_event_ind[unique_station_ind]][chan_iter]
                #Is it East?
                if channel[2]=='E':
                    chan_E=predparam_r[unique_event_ind[unique_station_ind]][chan_iter]
                    E_counter=E_counter+1
                elif channel[2]=='N':
                    N_counter=N_counter+1
                    chan_N=predparam_r[unique_event_ind[unique_station_ind]][chan_iter]
            
            #If for this recording (event and station combo) there is both an E and N
            #   reading, compute the geometrical average:
            if E_counter+N_counter==2:
                pga_recording_i=10**((log10(chan_E)+log10(chan_N))/2)

                #Now that the geometric avreage has been found for this reording, append
                #   the recording to the lists:
                #Save the event nuber as an integer:
                evnum.append(int(evnum_r[unique_event_ind[unique_station_ind]][chan_iter]))
                evlat.append(evlat_r[unique_event_ind[unique_station_ind]][chan_iter])
                evlon.append(evlon_r[unique_event_ind[unique_station_ind]][chan_iter])
                evdep.append(evdep_r[unique_event_ind[unique_station_ind]][chan_iter])
                sta.append(sta_r[unique_event_ind[unique_station_ind]][chan_iter])
                stlat.append(stlat_r[unique_event_ind[unique_station_ind]][chan_iter])
                stlon.append(stlon_r[unique_event_ind[unique_station_ind]][chan_iter])
                stelv.append(stelv_r[unique_event_ind[unique_station_ind]][chan_iter])
                grcircle.append(grcircle_r[unique_event_ind[unique_station_ind]][chan_iter])
                ml.append(ml_r[unique_event_ind[unique_station_ind]][chan_iter])
                mw.append(mw_r[unique_event_ind[unique_station_ind]][chan_iter])
                predparam_snr.append(predparam_snr_r[unique_event_ind[unique_station_ind]][chan_iter])
                
                if get_year=='yes':
                    evdate.append(evdate_r[unique_event_ind[unique_station_ind]][chan_iter])
                
                #Save the pga as the current geometical average:
                predparam.append(pga_recording_i)
    
    #At the end, convert the lists to arrays:
    evnum=array(evnum)
    evlat=array(evlat)
    evlon=array(evlon)
    evdep=array(evdep)
    sta=array(sta)
    stlat=array(stlat)
    stlon=array(stlon)
    stelv=array(stelv)
    grcircle=array(grcircle)
    ml=array(ml)
    mw=array(mw)
    predparam=array(predparam)
    predparam_snr=array(predparam_snr)
    
    if get_year=='yes':
        evdate=array(evdate)
        # split to get the year out of the date:
        evyear = zeros(len(evdate))
        for recordingi in range(len(evdate)):
            yeari = evdate[recordingi].split('/')[2]
            evyear[recordingi]=yeari
    
    ###Get indices for event and station for the sources.in and receivers.in files####
    ##Events first:
    
    #Get the unique station and event indices:
    unique_events=unique(evnum)
    
    #Zero out source ind array:
    source_ind=zeros((len(evnum)))
    #For each event in the record, devent, give it the source index to be used:
    for event_ind in range(len(unique_events)):
        eventi=unique_events[event_ind]
        
        #Find where in the recordings list the event number is the same as this one:
        recording_event_ind=where(evnum==eventi)
        
        #Set the source ind to be one plus this event, so it indexes with the raytracing program:
        source_ind[recording_event_ind]=event_ind+1
    
    #Now set these to integers...
    source_ind=source_ind.astype('int64')
    
    ##Next stations:
    unique_stations=unique(sta)
    
    #Zero out array:
    receiver_ind=zeros((len(sta)))
    #Loop through the unique stations:
    for station_ind in range(len(unique_stations)):
        stationi=unique_stations[station_ind]
        
        #Find where in the recordings list the station is the same as this one:
        recording_station_ind=where(sta==stationi)[0]
    
        #Set the receiver ind to be one plus this station, so it indexes with the raytracin gprogram:
        receiver_ind[recording_station_ind]=station_ind+1    
        
    #Set these to integers:
    receiver_ind=receiver_ind.astype('int64')
    
    ######
    #At the end, convert the lists to arrays:
    source_i=source_ind
    receiver_i=receiver_ind
    
    # If year was requested:
    if get_year=='yes':
        #Return the data - if PGA, predparam is pga in mgal, if PGV, predparam is pgv in cm/s:
        return evnum,evlat,evlon,evdep,sta,stlat,stlon,stelv,grcircle,ml,mw,predparam,source_i,receiver_i,predparam_snr,evyear
        
    elif get_year=='no':
        #Return the data - if PGA, predparam is pga in mgal, if PGV, predparam is pgv in cm/s:
        return evnum,evlat,evlon,evdep,sta,stlat,stlon,stelv,grcircle,ml,mw,predparam,source_i,receiver_i,predparam_snr

    
######
#Get Rrup
def compute_rrup(evlon,evlat,evdepth,stlon,stlat,stelv):
    '''
    Compute Rrup given the event and station lon,lat,z positions - ONLY USE ON POINT SOURCES!
    Input:
        evlon:          Array with event longitudes (deg)
        evlat:          Array with event latitudes (deg)
        evdepth:        Array with event depths (km)
        stlon:          Array with station longitudes (deg)
        stlat:          Array with station latitudes (deg)
        stdepth:        Array with station depths (km)
    Output:
        Rrup:           Array with Rrup distances (km)
    '''
    
    from pyproj import Proj
    from numpy import sqrt
    
    #Convert the event and station latitude/longitude to UTM x and y:
    #Make the projection:
    p=Proj(proj='utm',zone='11S',ellps='WGS84',inverse=True)
    
    #Convert the latitude and longitude to UTM X and Y (in meters)    
    evx,evy=p(evlon,evlat)
    stx,sty=p(stlon,stlat)
    
    #Event and station depth is currently in km; convert to m:
    evz=evdepth*1000
    #stations are negative, have positive depth:
    stz=stelv*-1000
    
    #Get distance Rrup - closest distance to rupture.
    #   Since they are almost all point sources, this can just be site/event distance:
    Rrup=sqrt((evx-stx)**2 + (evy-sty)**2 + (evz-stz)**2)    
    
    #Convert back to km:
    Rrup=Rrup/1000
    
    #Return:
    return Rrup
    
######
#Interpolate California Vs30 model for stations
def interp_vs30(stlat,stlon,vs30ascii):
    '''
    Interpolate the vs30 ascii file for vs30 values at select stations
    Input:
        sta:        List or array with strings of station names
        stlat:      Array with station latitudes
        stlon:      Array with station longitudes
        vs30ascii:  String with path to the vs30 model, with no header and columns:
                        long  lat  vs30
    Output:
        vs30:       Array with vs30 values at those stations
    '''
    
    from numpy import genfromtxt,meshgrid,unique
    from numpy import genfromtxt,sqrt,zeros,argmin
       
    #Import the ascii vs30 model:
    vs30_dat=genfromtxt(vs30ascii)
    x=vs30_dat[:,0]
    y=vs30_dat[:,1]
    z=vs30_dat[:,2]
    
    #Make the x and y meshgrid:
    X,Y=meshgrid(unique(x),unique(y))
    
    #Set vs30 array:
    vs30=zeros(len(stlon))
    
    #Get minimum distance (in degrees) - for every station, find the distance to
    #    every point in the model; find the minimum distance, adn this is where
    #    to take the vs30 value:
    for stai in range(len(stlon)):
        dist=sqrt((stlon[stai]-x)**2 + (stlat[stai]-y)**2)
        vs30[stai]=z[argmin(dist)]
            
    #Return:
    return vs30


########
# Compare Proxy Vs30 ID from Alan Yong's R script to the value from his 2016 paper
def vs30proxy_id2vs30(vs30_idfile,vs30_conversionfile,vs30_outfile):
    '''
    Input:
        vs30_idfile:                Path to the vs30 ID file from Yong's R script (sta, lon, lat, vs30 id)
        vs30_conversionfile:        Path to the vs30 conversion code file
        vs30_outfile:               Path to the outputfile for vs30
    Output: 
        vs30_outfile:               File with station name, longitude, latitude, proxy-based Vs30
    '''
    
    import numpy as np
    from string import replace
        
    # First read in the idfile:
    id_sta = np.genfromtxt(vs30_idfile,skip_header=1,usecols=1,dtype='S')
    id_data = np.genfromtxt(vs30_idfile,skip_header=1,usecols=range(2,5))
    
    # Get rid of the " in id_
    for k in range(len(id_sta)):
        tmpstring = replace(id_sta[k],'"','')
        id_sta[k]=tmpstring
    
    # Also read in the conversion code file:
    conversion_code = np.genfromtxt(vs30_conversionfile,skip_header=1)
    
    # Start a new vector, vs30 proxy, to save with id_sta and id_data:
    out_proxyVs30 = np.zeros(len(id_sta))
    
    # for every id in id_data, find the vs30 that corresponds to it:
    for sitei in range(len(id_sta)):
        # What is the ID to search for?
        id_i = id_data[sitei,2]
        conversion_ind = np.where(conversion_code[:,0]==id_i)[0][0]
        
        # Get the Vs30 that corresponds to that index:
        out_proxyVs30[sitei] = conversion_code[conversion_ind,1]
    
    # Now save the data to a file, write line by line:
    out_header = 'Sta \t Lon \t Lat \t Vs30 \n'
    
    outfile = open(vs30_outfile,'w')
    
    outfile.write(out_header)
    
    for linei in range(len(id_sta)):
        line_out = '%s\t%12.8f\t%10.8f\t%5.1f\n' % (id_sta[linei],id_data[linei,0],id_data[linei,1],out_proxyVs30[linei])
        outfile.write(line_out)
        
    outfile.close()
    
    
##############
def get_measured_vs30(vs30_idfile,measured_vs30file,output_measuredVs30file):
    '''
    Input:
        vs30_idfile:                Path to the idfile with at least the first three columns: sta, lon, lat - AND NO HEADER!!
        measured_vs30file:          Path to the CSV file from https://earthquake.usgs.gov/data/vs30/us/
        output_measuredVs30file:    Path to the output file with measured Vs30 values
    Output:
        output_measuredVs30file:    File with the columns: Sta, lon, lat, 
    '''
    
    import numpy as np
    from string import replace

    
    # Get the first column as a string, this has who measured it:
    surveyor_site = np.genfromtxt(measured_vs30file,skip_header=1,usecols=0,dtype='S',delimiter=',')
    
    surveyor = []
    site = []
    
    for elementi in range(len(surveyor_site)):
        surveyor.append(surveyor_site[elementi].split('.')[0])
        site.append(surveyor_site[elementi].split('.')[1])
    
    
    # Also get the method used:
    vs30_method = np.genfromtxt(measured_vs30file,skip_header=1,usecols=6,dtype='S',delimiter=',')
    # And the Vs30:
    vs30_measured = np.genfromtxt(measured_vs30file,skip_header=1,usecols=7,delimiter=',')
    
    
    # Get the station names/lon lat to use for printing:
    sta = np.genfromtxt(vs30_idfile,skip_header=1,usecols=0,dtype='S')
    sta_dat = np.genfromtxt(vs30_idfile,skip_header=1,usecols=range(1,3))
        
    # Get rid of the " in sta:
    for k in range(len(sta)):
        replace(sta[k],'"','')
    
    # Now for every site in the database, see if there is a measured value of Vs30 for it.
    # first make a vector to save them in:
    sta_out = []
    lon_out = []
    lat_out = []
    vs30_measured_out = []
    vs30_method_out = []
    
    # then loop through the stations provided to see if there's a measured value:
    for station in range(len(sta)):
        stationi = sta[station]
        
        measured_ind = np.where(np.array(site)==stationi)[0]
        
        # If an entry exists:
        if len(measured_ind)>0:
            # And if Alan Yong measured it, then grab the data:
            if (surveyor[measured_ind]=='AY'):
                sta_out.append(sta[station])
                lon_out.append(sta_dat[station,0])
                lat_out.append(sta_dat[station,1])
                vs30_measured_out.append(vs30_measured[measured_ind])
                vs30_method_out.append(vs30_method[measured_ind])
            
    
    ## Write out:
    outfile = open(output_measuredVs30file,'w')
    outfile.write('Sta \t Lon \t Lat \t Vs30 \t Method \n')
    
    for station in range(len(sta_out)):
        outfile.write('%s\t%12.8f\t%10.8f\t%5.1f\t%s\n' % (sta_out[station], lon_out[station], lat_out[station], vs30_measured_out[station], vs30_method_out[station][0]))
    
    outfile.close()
    
#############
# Combine measured and proxy stations for a given list of stations
def combine_measured_proxy_vs30(vs30_proxyfile,vs30_measuredfile,vs30_combinedfile):
    '''
    Combine measured and proxy Vs30's into one file with a flag.  Prioritize measured vs30.
    Input:
        vs30_proxyfile:             Path to the file with proxy vs30
        vs30_measuredfile:          Path to the file with measured Vs30
        vs30_combinedfile:          Path to the file with combined VS30
    Output:
        vs30_combinedfile:          File with output combined Vs30:  sta  lon  lat  vs30  method_flag(method, or proxy)
        station:                    Array with strings of station names
        vs30_out:                   Array with the Vs30
        vs30_method_out:            Array with a string of the method type - proxy, or the specified method
    '''
    
    import numpy as np
    from string import replace
    
    # The proxy file has all of the stations (or should), so use these as a starting point, read them in as the station lat and lon list:
    sta = np.genfromtxt(vs30_proxyfile,skip_header=1,usecols=0,dtype='S')
    sta_data = np.genfromtxt(vs30_proxyfile, skip_header=1,usecols=range(1,4))
    
    # Get rid of the " in sta:
    for k in range(len(sta)):
        replace(sta[k],'"','')
    
    # Read in also the measured info:
    sta_measured = np.genfromtxt(vs30_measuredfile,skip_header=1,usecols=0,dtype='S')
    sta_measured_vs30 = np.genfromtxt(vs30_measuredfile,skip_header=1,usecols=3)
    sta_measured_method = np.genfromtxt(vs30_measuredfile,skip_header=1,usecols=4,dtype='S')
        
    # Get rid of the " in sta:
    for k in range(len(sta_measured)):
        replace(sta_measured[k],'"','')
          
    # Set up empty arrays for all these things to go in:
    sta_out=[]
    lon_out=np.zeros(len(sta))
    lat_out=np.zeros(len(sta))
    vs30_out=np.zeros(len(sta))
    vs30_method_out=[]
    
    # For every station, first check if there's a measured value:
    for stationi in range(len(sta)):
        # if there's a measured value for this station, what's the index in the measured file:
        measured_ind = np.where(sta_measured==sta[stationi])[0]
        
        # Set the station and lon/lat:
        sta_out.append(sta[stationi])
        lon_out[stationi] = sta_data[stationi,0]
        lat_out[stationi] = sta_data[stationi,1]
        
        if len(measured_ind)>0:
            vs30_out[stationi] = sta_measured_vs30[measured_ind]
            vs30_method_out.append(sta_measured_method[measured_ind][0])
        else:
            vs30_out[stationi] = sta_data[stationi,2]
            vs30_method_out.append('proxy')
    
    # Convert sta to an array:
    sta_out = np.array(sta_out)
    vs30_method_out = np.array(vs30_method_out)
        
    # Write them to a file:
    outfile = open(vs30_combinedfile,'w')
    #writ header:
    outfile.write('Sta\tLon\tLat\tVs30\tMethod\n')
    
    for k in range(len(sta_out)):
        outfile.write('%s\t%12.8f\t%10.8f\t%5.1f\t%s\n' % (sta_out[k],lon_out[k],lat_out[k],vs30_out[k],vs30_method_out[k]))
    outfile.close()
    
    # Also return:
    return sta_out,vs30_out,vs30_method_out
    
    
############################################################################
def match_pga_pgv(evnum_pga,evlat,evlon,evdep,sta_pga,stlat,stlon,stelv,grcircle,ml,mw,pga_millig,pga_snr,evnum_pgv,sta_pgv,pgv_cmsec,pgv_snr,evyear='None'):
    '''
    Match PGA and PGV databases to have the same recordings, and
    be sorted in the same order.
    Input:
        evnum_pga:              Array with event numbers for PGA database
        evlat:                  Array with event lats for PGA database
        evlon:                  Array with event lons for PGA database
        evdep:                  Array with event depths for PGA database
        sta_pga:                Array with station names for PGA database
        stlat:                  Array with station lats for PGA database
        stlon:                  Array with station lons for PGA database
        stelv:                  Array with station elevs for PGA database
        grcircle:               Array with great circle path for PGA database
        ml:                     Array with local mag for PGA database
        mw:                     Array with moment mag for PGA database
        pga_millig:             Array with PGA in millig for PGA database
        pga_snr:                Array with PGA signal to noise for PGA database
        evnum_pgv:              Array with event numbers for PGV database
        sta_pgv:                Array with station names for PGV database
        pgv_cmsec:              Array with PGV in cm/sec for PGV database
        pgv_snr:                Array with PGV signal to noise for PGV database
        evyear:                 If not 'None', array with year of origin time for event
    Output:
    '''
    
    import numpy as np
    
    # Concatenate event numbers and station names into strings, to use as a unique identifier:
    evnum_pga_str = evnum_pga.astype('str')
    station_pga_str = sta_pga.astype('str')
    
    evnum_pgv_str = evnum_pgv.astype('str')
    station_pgv_str = sta_pgv.astype('str')
    
    # Add together element wise:
    recording_pga_str = np.core.defchararray.add(evnum_pga_str,station_pga_str)
    recording_pgv_str = np.core.defchararray.add(evnum_pgv_str,station_pgv_str)
    
    # First check that pga and pgv are unique sets:
    if len(np.unique(recording_pga_str))==len(recording_pga_str):
        print('PGA is a unique set of recordings')
    else:
        print('WARNING!!!! PGA IS NOT UNIQUE!!!')
        
    if len(np.unique(recording_pgv_str))==len(recording_pgv_str):
        print('PGV is a unique set of recordings')
    else:
        print('WARNING!!!! PGV IS NOT UNIQUE!!!')
    
    # Find the boolean intersection - first for where in pga pgv lies, then the opposite:
    intersect_array_pga = np.in1d(recording_pga_str,recording_pgv_str)
    intersect_array_pgv = np.in1d(recording_pgv_str,recording_pga_str)
    
    # Get the indices where they match - for pga, then pgv:
    match_indices_pga = np.where(intersect_array_pga==True)[0]
    match_indices_pgv = np.where(intersect_array_pgv==True)[0]
    
    # Check that the matched sets match in event number, line by line.
    # Set a counter to 0 before entering - if it's changed to 2, can continue; else won't:
    checkmatch_counter = 0
    if len(np.where((evnum_pga[match_indices_pga] - evnum_pgv[match_indices_pgv])==0)[0])==len(evnum_pga[match_indices_pga]):
        print('Events match...')
        checkmatch_counter+=1
    else:
        print('Events do not match, check something...')
        
    if len(set(sta_pga[match_indices_pga]) - set(sta_pgv[match_indices_pgv]))==0:
        print('Stations match...')
        checkmatch_counter+=1
    else:
        print('Stations do not match, check something...')
    
    # If they both match, checkmatch_counter will be 2:
    if checkmatch_counter==2:
        
        # these are all from the PGA database, so will use match_indices_pga:
        evnum = evnum_pga[match_indices_pga]
        evlat = evlat[match_indices_pga]
        evlon = evlon[match_indices_pga]
        evdep = evdep[match_indices_pga]
        sta = sta_pga[match_indices_pga]
        stlat = stlat[match_indices_pga]
        stlon = stlon[match_indices_pga]
        stelv = stelv[match_indices_pga]
        grcircle = grcircle[match_indices_pga]
        ml = ml[match_indices_pga]
        mw = mw[match_indices_pga]
        pga_millig = pga_millig[match_indices_pga]
        pga_snr = pga_snr[match_indices_pga]
        
        # If evyear is given (not None), then grab it here:
        if evyear!='None':
            evyear = evyear[match_indices_pga]
        
        # these are from PGV databse, so use match_indices_pgv:
        pgv_cmsec = pgv_cmsec[match_indices_pgv]
        pgv_snr = pgv_snr[match_indices_pgv]
        
        ####
        ## Now will need to get new source_i and receiver_i:
        source_i,receiver_i = find_source_receiver_indices(evnum,sta)
        
        # Then return:
        if evyear=='None':
            return evnum,evlat,evlon,evdep,sta,stlat,stlon,stelv,grcircle,ml,mw,pga_millig,pga_snr,pgv_cmsec,pgv_snr,source_i,receiver_i
        else:
            return evnum,evlat,evlon,evdep,sta,stlat,stlon,stelv,grcircle,ml,mw,pga_millig,pga_snr,pgv_cmsec,pgv_snr,source_i,receiver_i,evyear